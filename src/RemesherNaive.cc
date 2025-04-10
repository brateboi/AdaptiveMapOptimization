
#include "RemesherNaive.hh"




namespace AdaptiveMapOptimizer {


/**
 * Concrete implementation of where to split and collapse edges
 *
*/
void RemesherNaive::split_needles_and_collapse_spikes(
    int& split_count,
    int& collapse_count,
    bool splits_allowed,
    bool collapses_allowed){
    std::cout<<" ===== Remeshing mesh by splitting all edges once "<<std::endl;


    EdgePriorityQueue split_queue;
    EdgePriorityQueue collapse_queue;



    int spike_needle_elements_count = 0;
    int unhandled_count = 0;


    for(const auto& edge: target_mesh.edges()){


        if(target_mesh.status(edge).deleted()){
            continue;
        }

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

                //std::cout<<" -> added edge "<<edge.v0()<<", "<<edge.v1()<<" with op angle "<<max_angle<<" to split queue"<<std::endl;
                return true;
            }
            return false;
        };


        bool split_candidate = check_split(dir1, dir2, dir3, dir4);
        if (split_candidate){
            unhandled_count++;
        }

    }


    std::cout<<" -> found "<<collapse_queue.size()<<" edges to collapse and "<<split_queue.size()<<" to split"<<std::endl;
    std::cout<<" -> of a total of "<<target_mesh.n_edges()<<" edges: "<<std::endl;
    std::cout<<"      spikes: "<<collapse_queue.size()<<std::endl;
    std::cout<<"      others: "<<unhandled_count<<std::endl;


    if (splits_allowed){
        split_count    = split_edges_based_on_priority(split_queue);
    }


    std::cout << "===== Done remeshing:" << std::endl;
    std::cout << "    Edge splits: " << split_count << std::endl;

    return;
}

} // end namespace AdaptiveMapOptimizer
