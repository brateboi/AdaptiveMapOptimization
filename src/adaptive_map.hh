#ifndef ADAPTIVE_MAP_HH
#define ADAPTIVE_MAP_HH
#include "RemesherNaive.hh"
#include "RemesherValentin.hh"
#endif // ADAPTIVE_MAP_HH

#include "includes.hh"



namespace AdaptiveMapOptimizer {

using OptimizationTarget = std::vector<std::pair<OM::VertexHandle, OM::Vec3d>>;



class AdaptiveMapOptimizer {

public :


    AdaptiveMapOptimizer(TM reference_mesh_, OptimizationTarget &target_){
        reference_mesh = reference_mesh_;
        target_mesh = reference_mesh_;
        //target = target_;

        add_target_position(target_);

        // request status for both to allow for collapses
        reference_mesh.request_edge_status();
        reference_mesh.request_halfedge_status();
        reference_mesh.request_vertex_status();
        reference_mesh.request_face_status();

        target_mesh.request_edge_status();
        target_mesh.request_halfedge_status();
        target_mesh.request_vertex_status();
        target_mesh.request_face_status();

        constrain_non_original_vertices(reference_mesh);
        constrain_non_original_vertices(target_mesh);

        average_size = total_area(reference_mesh) / reference_mesh.n_faces();

    }


    // remeshes
    void remesh();

    void optimize_target_position();

    void one_step();

    void export_state(int state);

    void regularize_reference_mesh();

private:

    void add_target_position(OptimizationTarget &target_);
    void calculate_energy_of_target_mesh();

    // instance variables
    TM reference_mesh;
    TM target_mesh;
    //OptimizationTarget target;

    OM::VPropHandleT<OM::Vec3d> target_positions;
    OM::VPropHandleT<OM::Vec3d> remaining;
    OM::VPropHandleT<bool> has_target;

    RemesherValentin remesher = RemesherValentin(reference_mesh, target_mesh, remaining, target_positions, has_target);
    //RemesherNaive remesher = RemesherNaive(reference_mesh, target_mesh, remaining );

    int steps = 0;
    double average_size;

};

} // namespace AdaptiveMapOptimizer


void hello_world_adopt();


