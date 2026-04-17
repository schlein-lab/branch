// BRANCH v0.4 — `branch project` subcommand declaration.
//
// Three-layer reference projection for somatic variant detection.

#pragma once

namespace branch::cli {

/// Run the `branch project` subcommand.
///
/// @param argc  Argument count (excluding "branch" itself)
/// @param argv  Argument vector (argv[0] is "project")
/// @return Exit code (0 = success, 2 = usage error, 3+ = runtime error)
int run_project(int argc, char** argv);

}  // namespace branch::cli
