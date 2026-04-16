// BRANCH v0.2 — CLI entry point.
//
// Usage:
//   branch <subcommand> [args...]
//
// Subcommands:
//   analyze   paralog-aware CN inference from mosdepth regions.bed

#include <iostream>
#include <string>
#include <string_view>
#include <vector>

namespace branch::cli {

int run_analyze(int argc, char** argv);
int run_assemble(int argc, char** argv);

}  // namespace branch::cli

namespace {

void print_usage(std::ostream& os) {
    os << "branch 0.2.0 — somatic-mosaicism-aware assembler (in development)\n"
          "\nUsage:\n"
          "  branch <subcommand> [args...]\n"
          "\nSubcommands:\n"
          "  analyze    paralog-aware CN inference from mosdepth regions\n"
          "  assemble   reads -> minimizer overlap -> lossless graph (GFA-1.2)\n"
          "  version    print version and exit\n"
          "  help       print this help\n";
}

}  // namespace

int main(int argc, char** argv) {
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
    std::cerr << "Unknown subcommand: " << sub << "\n\n";
    print_usage(std::cerr);
    return 2;
}
