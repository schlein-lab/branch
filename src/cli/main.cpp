// BRANCH v0.2 — CLI entry point.
//
// Usage:
//   branch <subcommand> [args...]
//
// Subcommands:
//   analyze   paralog-aware CN inference from mosdepth regions.bed

#include <exception>
#include <iostream>
#include <new>
#include <string>
#include <string_view>
#include <vector>

#include "common/memory.hpp"

namespace branch::cli {

int run_analyze(int argc, char** argv);
int run_assemble(int argc, char** argv);
int run_project(int argc, char** argv);

}  // namespace branch::cli

namespace {

void print_usage(std::ostream& os) {
    os << "branch 0.2.0 — somatic-mosaicism-aware assembler (in development)\n"
          "\nUsage:\n"
          "  branch <subcommand> [args...]\n"
          "\nSubcommands:\n"
          "  analyze    paralog-aware CN inference from mosdepth regions\n"
          "  assemble   reads -> minimizer overlap -> lossless graph (GFA-1.2)\n"
          "  project    three-layer reference projection (v0.4 in development)\n"
          "  version    print version and exit\n"
          "  help       print this help\n";
}

}  // namespace

// Exit code 9 reserved for "memory budget exhausted" so shell wrappers
// and SLURM job scripts can distinguish OOM from other failures.
constexpr int kExitOutOfMemory = 9;

int dispatch(int argc, char** argv) {
    if (argc < 2) {
        print_usage(std::cout);
        return 0;
    }
    std::string_view sub = argv[1];
    if (sub == "help" || sub == "--help" || sub == "-h") {
        print_usage(std::cout);
        return 0;
    }
    if (sub == "version" || sub == "--version" || sub == "-V") {
        std::cout << "branch 0.2.0\n";
        return 0;
    }
    if (sub == "analyze") {
        return branch::cli::run_analyze(argc - 1, argv + 1);
    }
    if (sub == "assemble") {
        return branch::cli::run_assemble(argc - 1, argv + 1);
    }
    if (sub == "project") {
        return branch::cli::run_project(argc - 1, argv + 1);
    }
    std::cerr << "Unknown subcommand: " << sub << "\n\n";
    print_usage(std::cerr);
    return 2;
}

int main(int argc, char** argv) {
    branch::common::install_peak_rss_reporter();
    try {
        return dispatch(argc, argv);
    } catch (const std::bad_alloc&) {
        std::cerr << "\n[branch] FATAL: out of memory (std::bad_alloc). "
                     "Either raise --max-memory / SLURM --mem, or run on a "
                     "smaller input region. Peak RSS at abort = "
                  << branch::common::format_bytes(
                         branch::common::peak_rss_bytes())
                  << "\n";
        return kExitOutOfMemory;
    } catch (const std::exception& e) {
        std::cerr << "\n[branch] FATAL: " << e.what() << "\n";
        return 1;
    } catch (...) {
        std::cerr << "\n[branch] FATAL: unknown exception\n";
        return 1;
    }
}
