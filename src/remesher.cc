#include "includes.hh"
#include "helpers.hh"
#include "remesher.hh"
#include <queue>

/**
 * @brief triangle_flip_condition
 * @param mesh_
 * @param heh to be collapsed
 * @param _epsilon minimum sidelength of triangle
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

        // if sidelength of triangle is smaller than _epsilon return false

        double val0 = q[0] - p[0];
        double val1 = r[0] - p[0];
        double val2 = q[1] - p[1];
        double val3 = r[1] - p[1];

        if (val0 * val3 - val2 * val1 < 0) {
            std::cout << "       TRIANGLE FLIPPED" << std::endl;
            return false;
        }

        // only check if we are "CREATING" a small sidelength
        if (sidelength(p,q) < _epsilon || sidelength(p,r) < _epsilon){
            std::cout << "       SMALL SIDELENGTH" << std::endl;
            return false;
        }

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
 * @brief one_refine_step
 *
 * 1. Splits edges that are opposite to large angles 140 degrees = theta
 * 2. Collapse short edges that have ratio > 1/8 = epsilon of longest vs shortest edge
 * 3. both
 *
 * @return
 */
TM one_refine_step(TM mesh_, double max_angle_, double max_ratio_){
    int split_count = 0;
    int collapse_count = 0;
    split_needles_and_collapse_spikes(mesh_, max_angle_, max_ratio_, split_count, collapse_count);
    return mesh_;
}

void refine_and_optimize(TM mesh_){

    TM orig_mesh = mesh_;

    //optimize_inner_vertices(mesh_);

    TM after_refine = one_refine_step(orig_mesh);


    //TM after_refine_and_optimize = optimize_loose_vertices(after_refine);
    TM after_refine_and_optimize = optimize_inner_vertices(after_refine);



    int refine_optimize_times = 3;

    TM multiple_pass_mesh = mesh_;

    for (int i = 0; i < refine_optimize_times; ++i){
       multiple_pass_mesh = one_refine_step(multiple_pass_mesh);
       multiple_pass_mesh = optimize_inner_vertices(multiple_pass_mesh);
    }



    after_refine.garbage_collection();
    after_refine_and_optimize.garbage_collection();
    multiple_pass_mesh.garbage_collection();

    calculate_energy_of_mesh(after_refine);
    calculate_energy_of_mesh(after_refine_and_optimize);
    calculate_energy_of_mesh(multiple_pass_mesh);


    // OM::IO::write_mesh(orig_mesh, "../../../meshes/orig_mesh.om");
    OM::IO::write_mesh(after_refine, "../../../meshes/after_refine.om");
    OM::IO::write_mesh(after_refine_and_optimize, "../../../meshes/after_refine_and_optimize.om");
    OM::IO::write_mesh(multiple_pass_mesh, "../../../meshes/multiple_pass_mesh.om");

}

int split_edges_based_on_priority(TM &mesh_, EdgePriorityQueue split_queue){

    std::cout << "===== Splitting edges based on priority" << std::endl;
    std::cout << "      using incident triangles barycenter and edge endpoints within mesh, and midpoint on boundary of the mesh" << std::endl;

    int split_count = 0;
    int skipped_splits = 0;
    while (!split_queue.empty()) {

        auto top = split_queue.top();
        split_queue.pop();

        OM::EdgeHandle edge = top.first;
        double value = top.second;

        if (mesh_.status(edge).deleted()){
            skipped_splits++;
            //std::cout << "  -> skipping edge split "<< edge << "because already deleted" << std::endl;
            continue;
        }

        //std::cout << "Handling Edge: " << edge << ", Value: " << top.second << std::endl;


        //std::cout << "  -> Splitting edge " << mesh_.from_vertex_handle(heh0) << ", " << mesh_.to_vertex_handle(heh0) << std::endl;
        split_count++;

        mesh_.split(edge, mesh_.calc_edge_midpoint(edge));


    }
    //std::cout << "Skipped " << skipped_splits << " edge splits " << std::endl;

    return split_count;
}


int collapse_edges_based_on_priority(TM &mesh_, EdgePriorityQueue collapse_queue){

    std::cout << "===== Collapsing edges based on priority" << std::endl;

    int collapse_count = 0;
    int skipped_collapses = 0;
    while (!collapse_queue.empty()) {

        auto top = collapse_queue.top();
        collapse_queue.pop();

        OM::EdgeHandle edge = top.first;
        double value = top.second;

        if (mesh_.status(edge).deleted()){
            skipped_collapses++;
            //std::cout << "  -> skipping edge collapse "<< edge << "because already deleted" << std::endl;
            continue;
        }

        //std::cout << "Handling Edge: " << edge << ", Value: " << top.second << std::endl;

        auto collapsible = is_collapse_okay(mesh_, edge);
        if (collapsible.size() != 0){
            collapse_count++;
            mesh_.collapse(collapsible[0]);
        } else {
            std::cerr << " -> Could not collapse edge " << edge << std::endl;
            skipped_collapses++;
        }

    }

    //std::cout << "Skipped " << skipped_collapses << " edge splits " << std::endl;

    return collapse_count;
}



void split_needles_and_collapse_spikes(TM &mesh_,
                                      double split_angle_thresh,
                                      double collapse_side_len_ratio,
                                      int& split_count,
                                      int& collapse_count){

    std::cout<<" ===== Remeshing mesh by collapsing or splitting its edges "<<std::endl;
    std::cout<<" -       using needle angle threshold "<<split_angle_thresh<<std::endl;
    std::cout<<" - and spike collapse ratio threshold "<<collapse_side_len_ratio<<std::endl;


    EdgePriorityQueue split_queue;
    EdgePriorityQueue collapse_queue;


    OM::EPropHandleT<bool> split_prop;
    mesh_.add_property(split_prop, "split");
    mesh_.property(split_prop).set_persistent(true);

    OM::EPropHandleT<bool> collapse_prop;
    mesh_.add_property(collapse_prop, "collapse");
    mesh_.property(collapse_prop).set_persistent(true);

    OM::EPropHandleT<bool> both_prop;
    mesh_.add_property(both_prop, "both");
    mesh_.property(both_prop).set_persistent(true);




    int spike_needle_elements_count = 0;
    int unhandled_count = 0;


    for(const auto& edge: mesh_.edges()){


        if(mesh_.status(edge).deleted()){
            continue;
        }

        std::cout << " ----- checking edge " << edge << std::endl;

        mesh_.property(split_prop, edge) = false;
        mesh_.property(collapse_prop, edge) = false;
        mesh_.property(both_prop, edge) = false;

        auto heh0 = edge.h0();
        auto heh1 = edge.h1();

        // get points of incident triangles
        auto ov0 = mesh_.point(mesh_.to_vertex_handle(mesh_.next_halfedge_handle(heh0)));
        auto ov1 = mesh_.point(mesh_.to_vertex_handle(mesh_.next_halfedge_handle(heh1)));
        auto ev0 = mesh_.point(mesh_.to_vertex_handle(heh0));
        auto ev1 = mesh_.point(mesh_.to_vertex_handle(heh1));

        auto base_dir = ev0-ev1;

        auto dir1 = ev0 - ov0;
        auto dir2 = ev1 - ov0;

        auto dir3 = ev0 - ov1;
        auto dir4 = ev1 - ov1;

        auto check_split = [&](OM::Vec3d dir1, OM::Vec3d dir2, OM::Vec3d dir3, OM::Vec3d dir4) -> bool {


            auto angle = [](OM::Vec3d dir1, OM::Vec3d dir2 ) -> double {
                return 180.0 * std::acos(dir1.normalized().dot(dir2.normalized())) / M_PI;
            };

            double max_angle = std::max(angle(dir1, dir2), angle(dir3, dir4));

            if (max_angle > split_angle_thresh){
                split_queue.push({edge, -max_angle});
                mesh_.property(split_prop, edge) = true;

                //std::cout<<" -> added edge "<<edge.v0()<<", "<<edge.v1()<<" with op angle "<<max_angle<<" to split queue"<<std::endl;
                return true;
            }
            return false;
        };

        auto check_collapse = [&](OM::Vec3d dir1, OM::Vec3d dir2, OM::Vec3d dir3, OM::Vec3d dir4, OM::Vec3d base_dir) -> bool {
            double len1 = dir1.norm();
            double len2 = dir2.norm();
            double len3 = dir3.norm();
            double len4 = dir4.norm();
            double base_len = base_dir.norm();

            double max_len = std::max(std::max(len1, len2), std::max(len3, len4));

            double len_ratio = base_len / max_len;

            if (len_ratio < collapse_side_len_ratio){
                collapse_queue.push({edge, len_ratio});
                mesh_.property(collapse_prop, edge) = true;

                //std::cout<<" -> added edge "<<edge.v0()<<", "<<edge.v1()<<" with len ratio "<<len_ratio<<" to collapse queue"<<std::endl;
                return true;
            }
            return false;
        };

        bool collapse_candidate, split_candidate;
        collapse_candidate = check_collapse(dir1, dir2, dir3, dir4, base_dir);
        split_candidate = check_split(dir1, dir2, dir3, dir4);

        if (collapse_candidate && split_candidate){ // edge wants to be split and collapsed at the same time
            spike_needle_elements_count++;
            mesh_.property(both_prop, edge) = true;
        }

        if (!(collapse_candidate || split_candidate)){ // edge does not change
            unhandled_count++;
        }

    }


    std::cout<<" -> found "<<collapse_queue.size()<<" edges to collapse and "<<split_queue.size()<<" to split"<<std::endl;
    std::cout<<" -> of a total of "<<mesh_.n_edges()<<" edges: "<<std::endl;
    std::cout<<"      spikes: "<<collapse_queue.size()<<std::endl;
    std::cout<<"     needles: "<<split_queue.size()<<std::endl;
    std::cout<<"        both: "<<spike_needle_elements_count<<std::endl;
    std::cout<<"      others: "<<unhandled_count<<std::endl;

    //split_count    = split_edges_based_on_priority(mesh_, split_queue);
    collapse_count = collapse_edges_based_on_priority(mesh_, collapse_queue);

    std::cout << "===== Done remeshing:" << std::endl;
    std::cout << "    Edge splits: " << split_count << std::endl;
    std::cout << " Edge collapses: " << collapse_count << std::endl;

    return;
}

