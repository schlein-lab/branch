#pragma once
//
// bubble_finder — topological detection of bubbles in the unitig graph.
//
// A bubble is a pair (or set) of alternative paths between two anchor
// nodes that have multiple incoming/outgoing edges. We do NOT classify
// which side is h1 vs h2 here — that's phase_engine's job. We only
// emit the topological structure.

#include "types.hpp"
#include <vector>

namespace branch::wg::phaser {

struct BubbleFinderOpts {
    // Maximum length of an alternative path inside a bubble (in unitigs).
    // Anything longer is treated as "complex region" and not paired up.
    int max_alt_path_unitigs = 16;

    // Bubbles whose anchors are too far apart on the read backbone are
    // suspicious; flag them rather than phase them.
    int max_bubble_span_bp = 500'000;
};

class BubbleFinder {
public:
    // Annotates `unitigs` in-place by setting their `bubble` field, and
    // emits bubble structures into `out_bubbles`. The graph topology is
    // read from the next_nodes/prev_nodes fields of UnitigNode.
    void find(std::vector<UnitigNode>& unitigs,
              std::vector<Bubble>& out_bubbles,
              const BubbleFinderOpts& opts = {});
};

}  // namespace branch::wg::phaser
