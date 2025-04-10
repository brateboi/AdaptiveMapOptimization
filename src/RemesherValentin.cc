
#include "RemesherValentin.hh"




namespace AdaptiveMapOptimizer {


/**
 * Concrete implementation of where to split and collapse edges
 *
*/
void RemesherValentin::split_needles_and_collapse_spikes(
    int& split_count,
    int& collapse_count,
    bool splits_allowed,
    bool collapses_allowed){
    std::cout<<" ===== Remeshing mesh by collapsing or splitting its edges "<<std::endl;
    std::cout<<" -       using needle angle threshold "<<split_angle_thresh<<std::endl;
    std::cout<<" - and spike collapse ratio threshold "<<collapse_side_len_ratio<<std::endl;


    EdgePriorityQueue split_queue;
    EdgePriorityQueue collapse_queue;



    OM::EPropHandleT<bool> split_prop;
    target_mesh.add_property(split_prop, "split");
    target_mesh.property(split_prop).set_persistent(true);

    OM::EPropHandleT<bool> collapse_prop;
    target_mesh.add_property(collapse_prop, "collapse");
    target_mesh.property(collapse_prop).set_persistent(true);

    OM::EPropHandleT<bool> both_prop;
    target_mesh.add_property(both_prop, "both");
    target_mesh.property(both_prop).set_persistent(true);


    int spike_needle_elements_count = 0;
    int unhandled_count = 0;


    for(const auto& edge: target_mesh.edges()){


        if(target_mesh.status(edge).deleted()){
            continue;
        }

        //std::cout << " ----- checking edge " << edge << std::endl;

        target_mesh.property(split_prop, edge) = false;
        target_mesh.property(collapse_prop, edge) = false;
        target_mesh.property(both_prop, edge) = false;

        auto heh0 = edge.h0();
        auto heh1 = edge.h1();

        // get points of incident triangles
        auto ov0 = target_mesh.point(target_mesh.to_vertex_handle(target_mesh.next_halfedge_handle(heh0)));
        auto ov1 = target_mesh.point(target_mesh.to_vertex_handle(target_mesh.next_halfedge_handle(heh1)));
        auto ev0 = target_mesh.point(target_mesh.to_vertex_handle(heh0));
        auto ev1 = target_mesh.point(target_mesh.to_vertex_handle(heh1));

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
                target_mesh.property(split_prop, edge) = true;

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
                target_mesh.property(collapse_prop, edge) = true;

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
            target_mesh.property(both_prop, edge) = true;
        }

        if (!(collapse_candidate || split_candidate)){ // edge does not change
            unhandled_count++;
        }

    }


    std::cout<<" -> found "<<collapse_queue.size()<<" edges to collapse and "<<split_queue.size()<<" to split"<<std::endl;
    std::cout<<" -> of a total of "<<target_mesh.n_edges()<<" edges: "<<std::endl;
    std::cout<<"      spikes: "<<collapse_queue.size()<<std::endl;
    std::cout<<"     needles: "<<split_queue.size()<<std::endl;
    std::cout<<"        both: "<<spike_needle_elements_count<<std::endl;
    std::cout<<"      others: "<<unhandled_count<<std::endl;

    if (collapses_allowed){
        collapse_count = collapse_edges_based_on_priority(collapse_queue);
    }
    if (splits_allowed){
        split_count    = split_edges_based_on_priority(split_queue);
    }


    std::cout << "===== Done remeshing:" << std::endl;
    std::cout << "    Edge splits: " << split_count << std::endl;
    std::cout << " Edge collapses: " << collapse_count << std::endl;

    return;
}

} // end namespace AdaptiveMapOptimizer
