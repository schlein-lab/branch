# BRANCH Code Review: Classify Layer (Teil 2/4)

**Reviewer**: service-zyrkel  
**Datum**: 2025-01-21  
**Scope**: `src/classify/{features.hpp, cascade.hpp}`

Status audit: 2026-04-17 — markers added per finding.

## Files

| File | LOC | Zweck |
|------|-----|-------|
| features.hpp | 42 | 12 Feature-Enums, FeatureVector (48B float array), BubbleCandidate struct |
| cascade.hpp | 85 | CascadeConfig, 4-Stage classify_one() mit early-exit |

## Critical Issues

| Datei:Zeile | Severity | Problem | Status |
|-------------|----------|----------|--------|
| cascade.hpp:61 | MEDIUM | `stage_depth_ratio()` gibt raw DepthRatio zurück, aber Threshold ist 1.8 — Semantik unklar: ist 1.8 = 180% Coverage? Sollte normalisiert sein [0,1] | ⚠️ Partial: Threshold auf 2.0 angehoben (FIX 1, cascade.hpp L36), Semantik-Kommentar "Dup ~2x diploid" ergaenzt. Die raw-DepthRatio-Rueckgabe bleibt — keine [0,1]-Normalisierung, aber Confidence normalisiert (siehe naechste Zeile). |
| cascade.hpp:78-80 | HIGH | Stage 2 setzt confidence=1.0f hart statt p2 zu nutzen — Confidence ist nicht kalibriert | ✅ Resolved: FIX 2 berechnet `conf = min(1, p2/(threshold*1.5))` (cascade.hpp L100). |
| cascade.hpp:90-92 | MEDIUM | Stage 4 nutzt SegdupAnnotationFlag als Prior, aber ignoriert RepeatAnnotationFlag — asymmetrisch | ✅ Resolved: FIX 3 fuegt `stage_repeat_annotation()` hinzu, Stage 4 nimmt `max(p4_segdup, p4_repeat)` (cascade.hpp L112-114). |
| cascade.hpp:37-39 | LOW | CascadeConfig::min_bubble_length_bp/min_coverage/max_recursion_depth werden NICHT in classify_one() genutzt — Dead Config | ⚠️ Partial: `min_bubble_length_bp` wird jetzt als Early-Guard geprueft (cascade.hpp L86). `min_coverage` explizit entfernt (Kommentar L80-84: "coverage guard removed"). `max_recursion_depth` bleibt Dead-Config, weil Recursion in `classify_one()` nicht verdrahtet. |

## Classification Logic Feedback

**Positiv:**
- ✅ Early-exit Cascade-Pattern korrekt implementiert
- ✅ 4 Stages in sinnvoller Reihenfolge: Flanken → Depth → Span → Prior
- ✅ Feature-Enum stabil für ML-Input-Index
- ✅ BubbleCandidate hat genomische Koordinaten (chrom_id, start, end)

**Logik-Probleme:**
- ✅ Resolved — **Stage 2 Threshold-Semantik**: DepthRatio ≥ 1.8 → Duplication ist biologisch falsch. Threshold auf 2.0 korrigiert (FIX 1).
- ⚠️ Partial — **Keine Mixed-Klasse-Erkennung**: classify_one() gibt jetzt `BubbleClass::Mixed` zurueck wenn Flank+Depth widersprechen (FIX 5, cascade.hpp L119-126). Der rekursive SNP-Cluster-Dekompositionsmechanismus existiert in `mixed_decomposer.{hpp,cpp}` (getestet in `test_mixed_decomposer.cpp`), ist aber nicht aus `classify_one()` aufgerufen.
- ✅ Resolved — **Confidence nicht kalibriert**: Stage 2 jetzt normalisiert (FIX 2). Stage 1, 3, 4 geben feature-werte direkt zurueck, die per Design in [0,1] liegen (Jaccard, span-ratio, annotation-flag).

**Fehlende Discriminatoren:**
- ⚠️ Partial — ReadSpanCoverageIQR (Feature 3) wird im `feature_extractor.cpp` (L153-154) befuellt, aber in Cascade-Stages nicht genutzt.
- ⬜ Open — RefAlignScoreLeft/Right (Features 4,5) werden nicht genutzt.
- ⬜ Open — PathLengthDifference (Feature 6) wird nicht genutzt.
- ⚠️ Partial — GcContentDivergence (Feature 7) wird im `feature_extractor.cpp` (L157-167) befuellt, aber in Cascade-Stages nicht genutzt.
- ⬜ Open — AlleleFrequencyEstimate (Feature 10) wird nicht genutzt.
- ✅ Resolved — BubbleLengthBp (Feature 11) wird in `classify_one()` L85-87 als Early-Guard genutzt.

## TODOs

1. **MIXED_CLASS**: Implementiere rekursive Dekomposition für Mixed-Fälle mit den 5 Termination-Garantien aus dem Header-Kommentar — ⚠️ Partial: Mixed-Emission + `mixed_decomposer` vorhanden, nicht verdrahtet.
2. **DEPTH_THRESHOLD**: Ändern auf ~2.0 (oder ML-trained in v0.3) mit SegDup-Ausnahme — ✅ Resolved (2.0); SegDup-Ausnahme noch nicht implementiert.
3. **CALIBRATION**: Alle Stage-Confidences auf [0,1] normalisieren via Platt-Scaling o.ä. — ⚠️ Partial: Stage 2 kalibriert, andere Stages durch Feature-Design bereits in [0,1], kein Platt-Scaling.
4. **USE_ALL_FEATURES**: Features 3-7, 10-11 in Cascade einbauen — ⚠️ Partial: Features 3, 7, 11 extrahiert/genutzt; 4-6, 10 noch offen (fuer LightGBM v0.3).
5. **DEAD_CONFIG**: min_bubble_length_bp etc. tatsächlich in classify_one() prüfen — ⚠️ Partial: min_bubble_length_bp geprueft; min_coverage entfernt; max_recursion_depth noch Dead-Config.

---
*Review 2/4 abgeschlossen. Nächster Scope: backend/backend_vtable.hpp, gpu/overlap_kernel.cuh*
