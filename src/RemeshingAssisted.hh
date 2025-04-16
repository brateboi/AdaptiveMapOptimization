// RemesherBase.hh
#ifndef REMESHING_ASSISTED_HH
#define REMESHING_ASSISTED_HH

#include "includes.hh"
#include "helpers.hh"

namespace RemeshingAssisted {


// Define a custom comparator for the priority queue (min-heap)
struct EdgeCompare {
    bool operator()(const std::pair<OpenMesh::EdgeHandle, double>& a,
                    const std::pair<OpenMesh::EdgeHandle, double>& b) {
        return a.second > b.second;  // Min-heap: smaller double values have higher priority
    }
};

using EdgePriorityQueue = std::priority_queue<std::pair<OM::EdgeHandle, double>, std::vector<std::pair<OM::EdgeHandle, double>>, EdgeCompare>;


class RemeshingAssisted {
public:
    RemeshingAssisted(TM& reference_mesh_,
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

    void do_stuff() {
        return;

    }

protected:
    TM& reference_mesh;
    TM& target_mesh;
    OM::VPropHandleT<OM::Vec3d> &remaining;
    OM::VPropHandleT<OM::Vec3d> &target_positions;
    OM::VPropHandleT<bool> &has_target;





    int split_edges_based_on_priority(EdgePriorityQueue split_queue){

        return 0;
    };


    int collapse_edges_based_on_priority(EdgePriorityQueue collapse_queue){

        return 0;
    };



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
    std::vector<OM::HalfedgeHandle> is_collapse_okay(TM &mesh_, OM::EdgeHandle eh, double _epsilon = 1e-10){

        auto heh0 = mesh_.halfedge_handle(eh, 0);
        auto heh1 = mesh_.halfedge_handle(eh, 1);

        OM::VPropHandleT<bool> fixed_prop;
        if (!mesh_.get_property_handle(fixed_prop, "fixed")){
            std::cerr << "COULD NOT ACCESS FIXED PROPERTY, exiting" << std::endl;
            exit(-1);
        }


        auto check_all_cases = [&](OM::HalfedgeHandle &heh) -> bool {

            if (!mesh_.is_collapse_ok(heh)){
                std::cout << "-> LINK CONDITION NOT OKAY" << std::endl;
                return false;
            }


            auto from_v = mesh_.from_vertex_handle(heh); // from_vertex moves and gets deleted
            auto to_v = mesh_.to_vertex_handle(heh); // to_vertex doesnt move

            // Only non-original vertices can move
            if (mesh_.property(fixed_prop, from_v)){ // from_vertex not allowed to move because original vertex
                std::cout << "-> COLLAPSE ORIGINAL VERTEX NOT ALLOWED" << std::endl;
                return false;
            }

            // Collapse not along boundary not okay
            if (mesh_.is_boundary(from_v) && !mesh_.is_boundary(mesh_.edge_handle(heh))){
                std::cout << "-> COLLAPSE NOT ALONG BOUNDARY" << std::endl;
                return false;
            }

            // Geometry Check, check for inverted triangles around one-ring of halfedge.

            if (!triangle_flip_condition(mesh_, heh, _epsilon)){
                std::cout << "-> TRIANGLE FLIP CONDIITOIN NOT MET" << std::endl;
                return false;
            }

            // passed all checks, collapsing this halfedge is okay
            return true;
        };

        std::vector<OM::HalfedgeHandle> result = {};

        if (check_all_cases(heh0)) result.push_back(heh0);
        if (check_all_cases(heh1)) result.push_back(heh1);


        //passed all tests
        return result;
    }


    /**
 * @brief triangle_flip_condition
 * @param heh to be collapsed
 * @param _epsilon minimum area of triangle
 * @return true, if everything ok
 */
    bool triangle_flip_condition(TM &mesh_, OM::HalfedgeHandle &heh, double _epsilon){
        // simulate collapse of halfedge and check if any incident triangle flips or produces very thin triangle, i.e. side length < epsilon

        // hehs that need not be checked
        auto heh1 = mesh_.next_halfedge_handle(heh);
        auto heh2 = mesh_.prev_halfedge_handle(mesh_.opposite_halfedge_handle(heh));

        auto from_v = mesh_.from_vertex_handle(heh);
        auto to_v = mesh_.to_vertex_handle(heh);
        auto pos0 = mesh_.point(to_v);

        auto sidelength = [&](TM::Point p, TM::Point q) -> double {
            //std::cout << "checking point P " << p << "  and q  " << q << std::endl;


            double x2 = std::pow(p[0] - q[0], 2);
            double y2 = std::pow(p[1] - q[1], 2);


            double result = std::sqrt(x2+y2);
            if (result < _epsilon){
                std::cout << "---0SIDELENGTH" << result << std::endl;
            }
            return result;
        };


        // returns false if p, q, r are not oriented ccw, or if sidelength < _epsilon
        auto well_oriented = [&](TM::Point p, TM::Point q, TM::Point r) -> bool {

            // if area of triangle is smaller than _epsilon return false

            double val0 = q[0] - p[0];
            double val1 = r[0] - p[0];
            double val2 = q[1] - p[1];
            double val3 = r[1] - p[1];

            if (val0 * val3 - val2 * val1 < _epsilon) {
                std::cout << "       TOO SMALL Determinant" << std::endl;
                return false;
            }

            // only check if we are "CREATING" a small sidelength
            // if (sidelength(p,q) < _epsilon || sidelength(p,r) < _epsilon){
            //     std::cout << "       SMALL SIDELENGTH" << std::endl;
            //     return false;
            // }

            return true; // oriented ccw and sidelength > _epsilon
        };

        // Iterate over outgoing halfedges, collect remaining two points, check triangle orientation
        for (auto o_heh_it = mesh_.voh_begin(from_v); o_heh_it.is_valid(); ++o_heh_it){
            auto opp_heh = mesh_.next_halfedge_handle(*o_heh_it);

            // check if from face not needing to be checked
            if (opp_heh == heh1 || opp_heh == heh2){
                continue;
            }

            auto pos1 = mesh_.point(mesh_.from_vertex_handle(opp_heh));
            auto pos2 = mesh_.point(mesh_.to_vertex_handle(opp_heh));

            // triangle needs to have ccw orientation of pos0, pos1, pos2

            if (!well_oriented(pos0, pos1, pos2)){
                return false;
            }

        }

        // passed all checks
        return true;
    }

};

} // namespace AdaptiveMapOptimizer

#endif // REMESHING_ASSISTED_HH
