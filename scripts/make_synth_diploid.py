#!/usr/bin/env python3
"""Generate a synthetic diploid + HiFi-like reads with KNOWN read origin.

Produces ground truth for the dual-phaser neural voter and an end-to-end smoke
input for `branch wg`:

  <out>.hap1.fa        haplotype 1 (the reference backbone)
  <out>.hap2.fa        haplotype 2 (het SNVs + optional tandem-dup blocks)
  <out>.reads.fastq.gz reads sampled from both haps; the read NAME encodes the
                       true source: "<id> hap=1|2 contig=<c> pos=<p>"

Each read carries its true haplotype, so after running the de-novo + anchored
phasers we can label every disagreement by which track matched truth.

Deterministic given --seed; scale --genome-kb / --coverage / --contigs for large
training volumes.
"""
import argparse
import gzip
import random


def rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(rng, seq, snv_rate):
    """Return (mutated_seq, sorted het-SNV positions). Positions are shared
    between hap1 (backbone) and hap2 since no length change is introduced."""
    s = list(seq)
    het = []
    for i in range(len(s)):
        if rng.random() < snv_rate:
            s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
            het.append(i)
    return "".join(s), het


def add_tandem_dup(rng, seq, n_dups, dup_len):
    """Insert a few tandem duplications (CNV signal the BRANCH bubbles target)."""
    for _ in range(n_dups):
        if len(seq) < dup_len + 2:
            break
        start = rng.randint(0, len(seq) - dup_len)
        unit = seq[start:start + dup_len]
        seq = seq[:start + dup_len] + unit + seq[start + dup_len:]
    return seq


def sample_reads(rng, contigs_by_hap, het_by_contig, coverage, read_len,
                 err_rate, out_gz):
    import bisect
    n_written = 0
    with gzip.open(out_gz, "wt") as f:
        for hap, contigs in contigs_by_hap.items():
            for cname, seq in contigs:
                if len(seq) < read_len:
                    continue
                het = het_by_contig.get(cname, [])
                n_reads = max(1, (len(seq) * coverage) // read_len)
                for _ in range(n_reads):
                    pos = rng.randint(0, len(seq) - read_len)
                    r = list(seq[pos:pos + read_len])
                    for i in range(len(r)):
                        if rng.random() < err_rate:
                            r[i] = rng.choice("ACGT")
                    # het SNVs overlapped by this read → the read carries real
                    # haplotype signal only when this is > 0 (coords aligned
                    # between haps when no length-changing edits are present).
                    lo = bisect.bisect_left(het, pos)
                    hi = bisect.bisect_left(het, pos + read_len)
                    n_het = hi - lo
                    rid = f"r{n_written}"
                    f.write(f"@{rid} hap={hap} contig={cname} pos={pos} het={n_het}\n")
                    f.write("".join(r) + "\n+\n" + ("I" * read_len) + "\n")
                    n_written += 1
    return n_written


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--genome-kb", type=int, default=200, help="bp per contig (kb)")
    ap.add_argument("--contigs", type=int, default=2)
    ap.add_argument("--coverage", type=int, default=25)
    ap.add_argument("--read-len", type=int, default=10000)
    ap.add_argument("--snv-rate", type=float, default=0.001, help="het SNV rate")
    ap.add_argument("--err-rate", type=float, default=0.001, help="HiFi error rate")
    ap.add_argument("--dups", type=int, default=1, help="tandem dups per hap2 contig")
    ap.add_argument("--dup-len", type=int, default=2000)
    ap.add_argument("--ref-divergence", type=float, default=0.0,
                    help="extra SNV rate making the ANCHORED references differ "
                         "from the read-source haps (simulates mapping to a "
                         "reference panel rather than to self; 0 = anchor==truth)")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    clen = args.genome_kb * 1000

    hap1, hap2, het_by_contig = [], [], {}
    for c in range(args.contigs):
        cname = f"ctg{c}"
        backbone = rand_seq(rng, clen)
        h1 = backbone
        h2, het = mutate(rng, backbone, args.snv_rate)
        het_by_contig[cname] = het
        if args.dups > 0:
            # Note: tandem dups shift hap2 coordinates, so the het overlap
            # emitted in read names is only exact when --dups 0. The trainset
            # array runs with --dups 0 for clean het-region labelling.
            h2 = add_tandem_dup(rng, h2, args.dups, args.dup_len)
        hap1.append((cname, h1))
        hap2.append((cname, h2))

    def write_fa(path, contigs, hap):
        with open(path, "w") as f:
            for cname, seq in contigs:
                f.write(f">{cname}_hap{hap}\n")
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i + 80] + "\n")

    write_fa(args.out + ".hap1.fa", hap1, 1)
    write_fa(args.out + ".hap2.fa", hap2, 2)

    # Anchored references the phaser maps against. With --ref-divergence > 0
    # these differ from the read-source haps (realistic: mapping to a panel,
    # not to self), so the anchored track makes real errors and the voter must
    # learn feature-conditional trust instead of "always anchored".
    if args.ref_divergence > 0:
        aref1 = [(c, mutate(rng, s, args.ref_divergence)[0]) for c, s in hap1]
        aref2 = [(c, mutate(rng, s, args.ref_divergence)[0]) for c, s in hap2]
    else:
        aref1, aref2 = hap1, hap2
    write_fa(args.out + ".aref1.fa", aref1, 1)
    write_fa(args.out + ".aref2.fa", aref2, 2)
    n = sample_reads(rng, {1: hap1, 2: hap2}, het_by_contig,
                     args.coverage, args.read_len, args.err_rate,
                     args.out + ".reads.fastq.gz")
    print(f"wrote {args.out}.hap{{1,2}}.fa and {n} reads "
          f"({args.contigs}×{args.genome_kb}kb, {args.coverage}x)")


if __name__ == "__main__":
    main()
