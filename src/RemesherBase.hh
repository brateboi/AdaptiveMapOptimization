// RemesherBase.hh
#ifndef REMESHER_BASE_HH
#define REMESHER_BASE_HH

#include "includes.hh"
#include "helpers.hh"

namespace AdaptiveMapOptimizer {


// Define a custom comparator for the priority queue (min-heap)
struct EdgeCompare {
    bool operator()(const std::pair<OpenMesh::EdgeHandle, double>& a,
                    const std::pair<OpenMesh::EdgeHandle, double>& b) {
        return a.second > b.second;  // Min-heap: smaller double values have higher priority
    }
};

using EdgePriorityQueue = std::priority_queue<std::pair<OM::EdgeHandle, double>, std::vector<std::pair<OM::EdgeHandle, double>>, EdgeCompare>;



template<typename Impl>
class RemesherBase {
public:
    RemesherBase(TM& reference_mesh_,
                 TM& target_mesh_,
                 OM::VPropHandleT<OM::Vec3d> &remaining_,
                 OM::VPropHandleT<OM::Vec3d> &target_,
                 OM::VPropHandleT<bool> &has_target_)
        :
        reference_mesh(reference_mesh_),
        target_mesh(target_mesh_),
        remaining(remaining_),
        target_positions(target_),
        has_target(has_target_)

    {}

    void one_refine_step() {
        int split_count = 0, collapse_count = 0;

        // use implementation from derived class
        static_cast<Impl*>(this)->split_needles_and_collapse_spikes(
            split_count, collapse_count, true, true);
    }

protected:
    TM& reference_mesh;
    TM& target_mesh;
    OM::VPropHandleT<OM::Vec3d> &remaining;
    OM::VPropHandleT<OM::Vec3d> &target_positions;
    OM::VPropHandleT<bool> &has_target;

    /**
 *
 *  EdgeCollapseCheck
 *  1. is_collapse_ok from OpenMesh, Connectivity check
 *  2. Geometry okay:
 *      2.1 BoundaryEdge: - Both Original vertices          --> Not okay
 *                        - One original, one non-original  --> Non-original can be moved
 *                        - Both non-original               --> Both directions ok
 *      2.2 Geometry    :
 *          Check incident triangles for invertedness, use determinant with Epsilon on side length
 *          from remaining HalfEdges
 *
 *
*/


    int split_edges_based_on_priority(EdgePriorityQueue split_queue){

        std::cout << "===== Splitting edges based on priority" << std::endl;
        std::cout << "      using incident triangles barycenter and edge endpoints within mesh, and midpoint on boundary of the mesh" << std::endl;

        int split_count = 0;
        int skipped_splits = 0;
        while (!split_queue.empty()) {

            auto top = split_queue.top();
            split_queue.pop();

            OM::EdgeHandle edge = top.first;
            double value = top.second;

            if (target_mesh.status(edge).deleted()){
                skipped_splits++;
                //std::cout << "  -> skipping edge split "<< edge << "because already deleted" << std::endl;
                continue;
            }

            //std::cout << "Handling Edge: " << edge << ", Value: " << top.second << std::endl;


            //std::cout << "  -> Splitting edge " << mesh_.from_vertex_handle(heh0) << ", " << mesh_.to_vertex_handle(heh0) << std::endl;
            split_count++;



            // If edge lies on boundary, need to constrain the split vertex
            // Do this before actually splitting, because edge has different connectivitiy after split
            if (target_mesh.is_boundary(edge)){

                auto from_v = target_mesh.from_vertex_handle(target_mesh.halfedge_handle(edge));
                auto to_v = target_mesh.to_vertex_handle(target_mesh.halfedge_handle(edge));

                std::cout << "Edge that was on boundary was split, report targets" << std::endl;



                OM::Vec3d from = target_mesh.property(target_positions, from_v);
                OM::Vec3d to = target_mesh.property(target_positions, to_v);

                std::cout << "##### HAS TARGET FROM " << target_mesh.property(has_target, from_v) << std::endl;
                std::cout << "##### HAS TARGET TO   " << target_mesh.property(has_target, to_v) << std::endl;
                std::cout << "##### TARGET FROM " << from << std::endl;
                std::cout << "##### TARGET TO   " << to << std::endl;

                // split in both target and reference mesh
                auto split_vh_target = target_mesh.split(edge, target_mesh.calc_edge_midpoint(edge));
                auto split_vh_reference = reference_mesh.split(edge, reference_mesh.calc_edge_midpoint(edge));

                // TODO since the split is done at the midpoint, we can just average it. If this changes, change interpolation here
                target_mesh.property(target_positions, split_vh_target) = (from+to)/2.;
                target_mesh.property(has_target, split_vh_target) = true;
                target_mesh.property(remaining, split_vh_target) = OM::Vec3d(0.);

            } else {
                // split in both target and reference mesh
                auto split_vh_target = target_mesh.split(edge, target_mesh.calc_edge_midpoint(edge));
                auto split_vh_reference = reference_mesh.split(edge, reference_mesh.calc_edge_midpoint(edge));

                target_mesh.property(has_target, split_vh_target) = false;
                target_mesh.property(remaining, split_vh_target) = OM::Vec3d(0.);
            }


        }
        //std::cout << "Skipped " << skipped_splits << " edge splits " << std::endl;


        return split_count;

    };


    int collapse_edges_based_on_priority(EdgePriorityQueue collapse_queue){

        std::cout << "===== Collapsing edges based on priority" << std::endl;

        int collapse_count = 0;
        int skipped_collapses = 0;
        while (!collapse_queue.empty()) {

            auto top = collapse_queue.top();
            collapse_queue.pop();

            OM::EdgeHandle edge = top.first;
            double value = top.second;

            if (target_mesh.status(edge).deleted()){
                skipped_collapses++;
                //std::cout << "  -> skipping edge collapse "<< edge << "because already deleted" << std::endl;
                continue;
            }

            //std::cout << "Handling Edge: " << edge << ", Value: " << top.second << std::endl;

            auto collapsible_target = is_collapse_okay(target_mesh, edge);
            auto collapsible_reference = is_collapse_okay(reference_mesh, edge);

            // check if both vectors contain the same halfedge, if yes, collapse this one
            if (collapsible_target.size() > 0 && collapsible_reference.size() > 0 ){
                auto heh_reference = collapsible_reference[0];
                auto heh_target = collapsible_target[0];

                if (heh_reference == heh_target){
                    collapse_count++;
                    target_mesh.collapse(heh_target);
                    reference_mesh.collapse(heh_reference);
                } else {
                    std::cerr << " -> REFERENCE AND TARGET target edge did not match " << edge << std::endl;
                    skipped_collapses++;
                }
            } else {
                std::cerr << " -> Could not collapse edge " << edge << std::endl;
                skipped_collapses++;
            }

        }

        //std::cout << "Skipped " << skipped_collapses << " edge splits " << std::endl;

        return collapse_count;

    };

};

} // namespace AdaptiveMapOptimizer

#endif // REMESHER_BASE_HH
