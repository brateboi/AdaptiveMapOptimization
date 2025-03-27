#ifndef REMESHER_HH
#define REMESHER_HH
#endif // REMESHER_HH

#include "includes.hh"
#include "helpers.hh"

bool is_collapse_okay(TM &mesh, OM::HalfedgeHandle eh);

TM one_refine_step(TM mesh_, double max_angle_ = 130, double max_ratio_ = 1./8.);

void refine_and_optimize(TM mesh_);



// Define a custom comparator for the priority queue (min-heap)
struct EdgeCompare {
    bool operator()(const std::pair<OpenMesh::EdgeHandle, double>& a,
                    const std::pair<OpenMesh::EdgeHandle, double>& b) {
        return a.second > b.second;  // Min-heap: smaller double values have higher priority
    }
};

using EdgePriorityQueue = std::priority_queue<std::pair<OM::EdgeHandle, double>, std::vector<std::pair<OM::EdgeHandle, double>>, EdgeCompare>;


int split_edges_based_on_priority(TM &mesh_, EdgePriorityQueue split_queue);
int collapse_edges_based_on_priority(TM &mesh_, EdgePriorityQueue collapse_queue);
void split_needles_and_collapse_spikes(TM &mesh_,
                                       double split_angle_thresh,
                                       double collapse_side_len_ratio,
                                       int& split_count,
                                       int& collapse_count, bool splits_allowed, bool collapses_allowed);
