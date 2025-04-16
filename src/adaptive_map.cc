#include <format>
#include <includes.hh>
#include <adaptive_map.hh>

namespace AMO = AdaptiveMapOptimizer;
using OptimizationTarget = std::vector<std::pair<OM::VertexHandle, OM::Vec3d>>;

namespace AdaptiveMapOptimizer {


void AdaptiveMapOptimizer::remesh(){

    remesher.one_refine_step();

}


void AdaptiveMapOptimizer::optimize_target_position_with_hard_constraints(){

    using SDE_Elem = SDEB2D_PH_AD;

    std::cout<<" ---------------------------------------------------------------------------------------- "<<std::endl;
    std::cout<<" *** trying to deform mesh now with hard constraints" << std::endl;


    //low weight so the pull has the most influence and this only really serves as a barrier
    const double Dirichlet_w(1), Dirichlet_barrier_w(0), Dirichlet_barrier_min(1e1);
    const int max_inner_iters(100);
    const int n_vars(target_mesh.n_vertices() * 2);

    std::cout << "#######nr of n_vars " << n_vars << std::endl;

    COMISO::FiniteElementProblem fe_problem(n_vars);
    std::vector<COMISO::LinearConstraint> eq_constraints;


    COMISO::FiniteElementSet<SDE_Elem> fe_barrier("barrier elements");

    //add the barrier elements
    for (auto fh: target_mesh.faces()) {

        //std::cout<<" --------------------------------------------------- handling face "<<fh<<std::endl;

        double w = 1;
        double sizing_scale = 1.0; //compute_sizing_scale(ch);

        typename SDE_Elem::VecI vi_barrier;
        typename SDE_Elem::VecC vc_barrier;
        typename SDE_Elem::VecV vx_barrier;

        std::vector<TM::VertexHandle> verts = std::vector(target_mesh.fv_begin(fh), target_mesh.fv_end(fh));

        //auto verts = mesh_.get_halfface_vertices(mesh_.halfface_handle(fh,0));


        for (int local_vid(0); local_vid < 3; local_vid++) {
            const auto &vh = verts[local_vid];

            int global_id = vh.idx() * 2;
            for (int i(0); i < 2; i++) {
                vi_barrier[2 * local_vid + i] = global_id + i;
            }

            for (int i(0); i < 2; i++) {
                fe_problem.x()[global_id + i] = reference_mesh.point(vh)[i];
                vx_barrier[2 * local_vid + i] = fe_problem.x()[global_id + i];
            }
        }

        // use reference element from reference mesh
        auto p0 = reference_mesh.point(verts[0]);
        auto p1 = reference_mesh.point(verts[1]);
        auto p2 = reference_mesh.point(verts[2]);

        std::vector<Vec2d> P_ref(3);
        P_ref[0] << p0[0],p0[1];
        P_ref[1] << p1[0],p1[1];
        P_ref[2] << p2[0],p2[1];


        auto setup_res = SDE_Elem::compute_constants(Dirichlet_w, Dirichlet_barrier_w, Dirichlet_barrier_min,
                                                   P_ref,
                                                   vc_barrier);

        //std::cout<<" - constants vector: "<<vc_barrier.transpose()<<std::endl;

        if (setup_res) {
            std::cout << " ERROR: constants couldn't be computed for face: " <<fh<< std::endl;
            return;
        }

        fe_barrier.instances().add_element(vi_barrier,
                                           vc_barrier);
    }
    fe_problem.add_set(&fe_barrier);



    // constrain the fixed vertices, currently as a hard constraint, TODO implement as penalty term
    for (const auto vh : target_mesh.vertices()){

        if (target_mesh.property(has_target, vh)){

            auto pos = target_mesh.property(target_positions, vh);

            COMISO::LinearConstraint::SVectorNC coeffs0, coeffs1;
            coeffs0.resize(n_vars);
            coeffs1.resize(n_vars);

            int var_idx = vh.idx() * 2;
            coeffs0.coeffRef(var_idx + 0) = 1;
            coeffs1.coeffRef(var_idx + 1) = 1;


            eq_constraints.push_back(COMISO::LinearConstraint(coeffs0, -pos[0]));
            eq_constraints.push_back(COMISO::LinearConstraint(coeffs1, -pos[1]));
            std::cout << " - constrained vertex " << vh << " to " << pos << std::endl;
        }
    }


    std::vector<COMISO::NConstraintInterface*> constraints_pointers;
    constraints_pointers.reserve(constraints_pointers.size());
    for( auto& c : eq_constraints) {
        constraints_pointers.push_back(&c);
    }

    std::cout<<" energies before optimization: "<<std::endl;
    fe_problem.print_objectives();

    COMISO::NPTiming fet_problem(&fe_problem);

    COMISO::NewtonSolver ns( 1e-1, 1e-9, max_inner_iters);
    ns.solve_infeasible_start(&fet_problem, constraints_pointers);
    //solver_iters += ns.iters();
    bool converged = ns.converged();

    std::cout<<" energies after optimization: "<<std::endl;
    fe_problem.print_objectives();
    //copy problem data back to mesh

    for(auto v: target_mesh.vertices()){
        int global_id = 2 * v.idx();

        typename OM::Vec3d new_pos;
        for(int i = 0; i<2; i++){
            new_pos[i] = fe_problem.x()[global_id + i];
        }
        target_mesh.set_point(v, new_pos);
    }


    //return updated mesh positions

    std::cout<<" --------------------------------------------------------- "<<std::endl;

}

void AdaptiveMapOptimizer::optimize_target_position_with_penalty_terms(){

    using SDE_Elem = SDEB2D_PH_AD;

    std::cout<<" ---------------------------------------------------------------------------------------- "<<std::endl;
    std::cout<<" *** trying to deform mesh now with penalty term" << std::endl;


    //low weight so the pull has the most influence and this only really serves as a barrier
    const double Dirichlet_w(1), Dirichlet_barrier_w(0), Dirichlet_barrier_min(1e1);
    const int max_inner_iters(100);
    const int n_vars(target_mesh.n_vertices() * 2);

    const double target_penalty_weight = 1e10; // TODO, SCALE INVARIANT
    //const double target_penalty_weight = 1/(average_size*average_size);

    std::cout << "#######nr of n_vars " << n_vars << std::endl;

    COMISO::FiniteElementProblem fe_problem(n_vars);
    std::vector<COMISO::LinearConstraint> eq_constraints;


    COMISO::FiniteElementSet<SDE_Elem> fe_barrier("barrier elements");

    //add the barrier elements
    for (auto fh: target_mesh.faces()) {

        //std::cout<<" --------------------------------------------------- handling face "<<fh<<std::endl;

        double w = 1;
        double sizing_scale = 1.0; //compute_sizing_scale(ch);

        typename SDE_Elem::VecI vi_barrier;
        typename SDE_Elem::VecC vc_barrier;
        typename SDE_Elem::VecV vx_barrier;

        std::vector<TM::VertexHandle> verts = std::vector(target_mesh.fv_begin(fh), target_mesh.fv_end(fh));

        //auto verts = mesh_.get_halfface_vertices(mesh_.halfface_handle(fh,0));


        for (int local_vid(0); local_vid < 3; local_vid++) {
            const auto &vh = verts[local_vid];

            int global_id = vh.idx() * 2;
            for (int i(0); i < 2; i++) {
                vi_barrier[2 * local_vid + i] = global_id + i;
            }

            for (int i(0); i < 2; i++) {
                fe_problem.x()[global_id + i] = reference_mesh.point(vh)[i];
                vx_barrier[2 * local_vid + i] = fe_problem.x()[global_id + i];
            }
        }

        // use reference element from reference mesh
        auto p0 = reference_mesh.point(verts[0]);
        auto p1 = reference_mesh.point(verts[1]);
        auto p2 = reference_mesh.point(verts[2]);

        std::vector<Vec2d> P_ref(3);
        P_ref[0] << p0[0],p0[1];
        P_ref[1] << p1[0],p1[1];
        P_ref[2] << p2[0],p2[1];


        auto setup_res = SDE_Elem::compute_constants(Dirichlet_w, Dirichlet_barrier_w, Dirichlet_barrier_min,
                                                     P_ref,
                                                     vc_barrier);

        //std::cout<<" - constants vector: "<<vc_barrier.transpose()<<std::endl;

        if (setup_res) {
            std::cout << " ERROR: constants couldn't be computed for face: " <<fh<< std::endl;
            return;
        }

        fe_barrier.instances().add_element(vi_barrier,
                                           vc_barrier);
    }
    fe_problem.add_set(&fe_barrier);


    // ****************************************************
    // Add the penalty elements for vertices with target positions
    // ****************************************************
    COMISO::FiniteElementSet<QP2D_AD_PH> fe_penalty("penalty elements");

    for (const auto vh : target_mesh.vertices()) {
        if (target_mesh.property(has_target, vh)) {
            // Get the target position
            auto target_pos = target_mesh.property(target_positions, vh);

            // Global index in the optimization variable vector
            int global_id = vh.idx() * 2;

            // Prepare the constant vector for the penalty term:
            // [target_penalty_weight, target_x, target_y]
            QuadraticPenaltyElement2D::VecC consts;
            consts[0] = target_penalty_weight;
            consts[1] = target_pos[0];
            consts[2] = target_pos[1];

            // Define the index vector for this vertex
            typename QP2D_AD_PH::VecI indices;
            indices[0] = global_id;
            indices[1] = global_id + 1;

            // Add the penalty element to the finite element set
            fe_penalty.instances().add_element(indices, consts);

            std::cout << " - added penalty for vertex " << vh << " with target " << target_pos << std::endl;
        }
    }

    fe_problem.add_set(&fe_penalty);


    // ****************************************************
    // Optimization: solve the unconstrained minimization problem
    // (Note that the penalty terms penalize deviations from the targets)
    // ****************************************************

    std::cout<<" Energies before optimization: "<<std::endl;
    fe_problem.print_objectives();

    COMISO::NPTiming fet_problem(&fe_problem);

    COMISO::NewtonSolver ns( 1e-1, 1e-9, max_inner_iters);


    // Unconstrained problem now, use normal Newton solve
    ns.solve(&fet_problem);
    //solver_iters += ns.iters();
    bool converged = ns.converged();

    std::cout<<" Energies after optimization: "<<std::endl;
    fe_problem.print_objectives();
    //copy problem data back to mesh

    for(auto v: target_mesh.vertices()){
        int global_id = 2 * v.idx();

        typename OM::Vec3d new_pos;
        for(int i = 0; i<2; i++){
            new_pos[i] = fe_problem.x()[global_id + i];
        }
        target_mesh.set_point(v, new_pos);
    }

    //return updated mesh positions

    std::cout<<" --------------------------------------------------------- "<<std::endl;

}

void AdaptiveMapOptimizer::add_target_position(OptimizationTarget &target){
    std::cout << "Adding optimization target " << std::endl;

    std::cout << "nr of verts " << target_mesh.n_vertices() << std::endl;

    // only add if not already existing
    if (!target_mesh.get_property_handle(target_positions, "target_positions")){
        target_mesh.add_property(target_positions, "target_positions");
    }

    if (!target_mesh.get_property_handle(remaining, "remaining_dirs")){
        target_mesh.add_property(remaining, "remaining_dirs");
        target_mesh.property(remaining).set_persistent(true);
    }


    if (!target_mesh.get_property_handle(has_target, "has_target")){
        target_mesh.add_property(has_target, "has_target");
        target_mesh.property(has_target).set_persistent(true);
        // make default false
        for (const auto v : target_mesh.vertices()){
            target_mesh.property(has_target, v) = false;
        }
    }



    std::cout << "target vertices: ";
    for (auto &[vh, pos] : target){
        std::cout << vh << " ";
        target_mesh.property(target_positions, vh) = pos;
        target_mesh.property(has_target, vh) = true;

    }
    std::cout << std::endl;
}

void AdaptiveMapOptimizer::calculate_energy_of_target_mesh(){

    auto SDE_Element = [&](OM::FaceHandle face) -> double {
        std::vector<OM::VertexHandle> verts;

        for (auto v_it = target_mesh.fv_begin(face); v_it.is_valid(); ++v_it){
            verts.push_back(*v_it);
        }

        assert(verts.size() == 3);

        Vec2d v0 = Vec2d(target_mesh.point(verts[0]).data());
        Vec2d v1 = Vec2d(target_mesh.point(verts[1]).data());
        Vec2d v2 = Vec2d(target_mesh.point(verts[2]).data());

        // use reference element from reference mesh
        auto p0 = reference_mesh.point(verts[0]);
        auto p1 = reference_mesh.point(verts[1]);
        auto p2 = reference_mesh.point(verts[2]);

        std::vector<Vec2d> P_ref(3);
        P_ref[0] << p0[0],p0[1];
        P_ref[1] << p1[0],p1[1];
        P_ref[2] << p2[0],p2[1];

        // get ref edge  matrix
        Mat2d E_ref;

        E_ref.col(0) = P_ref[1] - P_ref[0];
        E_ref.col(1) = P_ref[2] - P_ref[0];


        // get edge vectors
        Mat2d E;
        E.col(0) = v1 - v0;
        E.col(1) = v2 - v0;


        // calc  Jacobi Matrix of map
        Mat2d J = E * E_ref.inverse();

        double E_sd = (J.squaredNorm() + J.inverse().squaredNorm() - 4.0);

        return E_sd;
        //if (barrier_w && barrier_min - E_sd <= 0.0) {
        //std::cout<<" inf barrier term"<<std::endl;
        //    return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
        //}
        //return w * E_sd;

        //return w * E_sd + barrier_w / (barrier_min - E_sd);

        //return 0;
    };

    // only add if not already existing
    OM::FPropHandleT<double> distortion;
    if (!target_mesh.get_property_handle(distortion, "e_distortion")){
        target_mesh.add_property(distortion, "e_distortion");
        target_mesh.property(distortion).set_persistent(true);
    } else {
        target_mesh.get_property_handle(distortion, "e_distortion");
    }

    for (const auto fh : target_mesh.faces()){
        target_mesh.property(distortion, fh) = SDE_Element(fh);
    }
}

void AdaptiveMapOptimizer::export_state(int state){

    target_mesh.garbage_collection();
    reference_mesh.garbage_collection();

    for (const auto vh : target_mesh.vertices()){
        if (target_mesh.property(has_target, vh)){
            auto pos = target_mesh.property(target_positions, vh);
            target_mesh.property(remaining, vh) = pos - target_mesh.point(vh);
        }
    }


    calculate_energy_of_target_mesh();


    std::string m1 = std::format("frames/{:02}target_mesh.om", state);
    std::string m2 = std::format("frames/{:02}reference_mesh.om", state);



    OM::IO::write_mesh(target_mesh, m1);
    OM::IO::write_mesh(reference_mesh, m2);
}

void AdaptiveMapOptimizer::regularize_reference_mesh(){

    using MyElem = SDEB2D_PH_AD;
    //using MyElem = COMISO::FiniteElementTinyAD<FoldoverFreeElement2D>;

    std::cout<<" --------------------------------------------------------- "<<std::endl;
    std::cout<<" *** regularizing reference mesh, only inner vertices of mesh scaling the reference element by average area"<<std::endl;


    std::vector<Vec2d> P_ref(3);
    P_ref[0] << 0,0;
    P_ref[1] << 1,0;
    P_ref[2] << 0.5, 0.5 * std::sqrt(3.0);


    // get average size of elements, 0.433013 area of reference element
    double scale_factor = std::sqrt(average_size / 0.433013);

    std::cout << "SCALE FACTOR " << scale_factor << std::endl;

    for(int i(0); i<3; i++){
        P_ref[i] *= scale_factor;
    }


    //low weight so the pull has the most influence and this only really serves as a barrier
    const double Dirichlet_w(1), Dirichlet_barrier_w(0), Dirichlet_barrier_min(1e8);
    const int max_inner_iters(100);
    const int n_vars(reference_mesh.n_vertices() * 2);

    COMISO::FiniteElementProblem fe_problem(n_vars);
    std::vector<COMISO::LinearConstraint> eq_constraints;

    COMISO::FiniteElementSet<MyElem> fe_barrier("barrier elements");

    //add the barrier elements
    for (auto fh: reference_mesh.faces()) {

        //std::cout<<" --------------------------------------------------- handling face "<<fh<<std::endl;

        double w = 1;
        double sizing_scale = 1.0; //compute_sizing_scale(ch);

        typename MyElem::VecI vi_barrier;
        typename MyElem::VecC vc_barrier;
        typename MyElem::VecV vx_barrier;

        //auto temp = mesh_.face(fh);
        std::vector<TM::VertexHandle> verts = std::vector(reference_mesh.fv_begin(fh), reference_mesh.fv_end(fh));

        //auto verts = mesh_.get_halfface_vertices(mesh_.halfface_handle(fh,0));


        for (int local_vid(0); local_vid < 3; local_vid++) {
            const auto &vh = verts[local_vid];

            int global_id = vh.idx() * 2;
            for (int i(0); i < 2; i++) {
                vi_barrier[2 * local_vid + i] = global_id + i;
            }

            for (int i(0); i < 2; i++) {
                fe_problem.x()[global_id + i] = reference_mesh.point(vh)[i];
                vx_barrier[2 * local_vid + i] = fe_problem.x()[global_id + i];
            }
        }

        // update P_ref each time to reference mesh

        auto setup_res = MyElem::compute_constants(Dirichlet_w, Dirichlet_barrier_w, Dirichlet_barrier_min,
                                                   P_ref,
                                                   vc_barrier);

        //std::cout<<" - constants vector: "<<vc_barrier.transpose()<<std::endl;

        if (setup_res) {
            std::cout << " ERROR: constants couldn't be computed for face: " <<fh<< std::endl;
            exit(-1);
        }

        fe_barrier.instances().add_element(vi_barrier,
                                           vc_barrier);
    }

    fe_problem.add_set(&fe_barrier);


    //constrain the fixed vertices
    for(auto v: reference_mesh.vertices()) {
        if(reference_mesh.is_boundary(v)) {
            COMISO::LinearConstraint::SVectorNC coeffs0, coeffs1;
            coeffs0.resize(n_vars);
            coeffs1.resize(n_vars);

            int var_idx = v.idx() * 2;
            coeffs0.coeffRef(var_idx + 0) = 1;
            coeffs1.coeffRef(var_idx + 1) = 1;

            auto pos = reference_mesh.point(v);


            eq_constraints.push_back(COMISO::LinearConstraint(coeffs0, -pos[0]));
            eq_constraints.push_back(COMISO::LinearConstraint(coeffs1, -pos[1]));
            //std::cout << " - constrained vertex " << v << " to " << pos << std::endl;
        }
    }

    std::vector<COMISO::NConstraintInterface*> constraints_pointers;
    constraints_pointers.reserve(constraints_pointers.size());
    for( auto& c : eq_constraints) {
        constraints_pointers.push_back(&c);
    }

    std::cout<<" energies before optimization: "<<std::endl;
    fe_problem.print_objectives();

    COMISO::NPTiming fet_problem(&fe_problem);


    COMISO::NewtonSolver ns( 1e-1, 1e-9, max_inner_iters);
    ns.solve_infeasible_start(&fet_problem, constraints_pointers);
    //solver_iters += ns.iters();
    bool converged = ns.converged();

    std::cout<<" energies after optimization: "<<std::endl;
    fe_problem.print_objectives();
    //copy problem data back to mesh

    for(auto v: reference_mesh.vertices()){
        int global_id = 2 * v.idx();

        typename OM::Vec3d new_pos;
        for(int i = 0; i<2; i++){
            new_pos[i] = fe_problem.x()[global_id + i];
        }
        reference_mesh.set_point(v, new_pos);
    }


    //return updated mesh positions


    std::cout<<" --------------------------------------------------------- "<<std::endl;

}

void AdaptiveMapOptimizer::one_step(){

    switch (constraint_mode) {
    case HardConstraint:

        optimize_target_position_with_hard_constraints();
        break;
    case ExactPenaltyTerm:

        std::cerr << "Not implemented error, exiting" << std::endl;
        exit(-1);

        break;

    case QuadraticPenaltyTerm:

        optimize_target_position_with_penalty_terms();
        break;

    default:
        std::cerr << "Undefined Constraint Mode, please specify" << std::endl;
        exit(-1);
        break;
    }

    export_state(steps++);
    remesh();
    regularize_reference_mesh();
    export_state(steps++);
}





void AdaptiveMapOptimizer::run_pipeline(){
    std::cout <<"Hello ???" << std::endl;
    switch (pipeline_mode){
    case Mine:


        one_step();
        one_step();
        one_step();
        one_step();
        std::cout << "lule" << std::endl;

        break;
    case RemeshingAssisted:



        std::cout << "Remeshing Assited" << std::endl;

        break;

    case Panozzo:
        std::cerr << "Not implemented error, exiting" << std::endl;
        exit(-1);

        break;
    default:
        std::cerr << "Undefined Pipeline Mode, please specify" << std::endl;
        exit(-1);

        break;

    }


    return;

}

} // namespace AdaptiveMapOptimizer




void hello_world_adopt()
{
    std::cout << "HELLO FROM MY LIBRARY" << std::endl;
}




int main(int argc, char const *argv[])
{
    std::cout << "Hello World " << std::endl;

    // TM mesh = get_disk_mesh_without_interior_vertices(10, 2);
    // OptimizationTarget target = get_disk_optimization_target(mesh, 10, 1);


    TM mesh = get_box_mesh_one_interior_vertex(30, 4);
    OptimizationTarget target = get_box_optimization_target(mesh, 4,2);

    //OM::IO::write_mesh(disk, "disk.om");
    //OM::IO::write_mesh(box, "box.om");

    //auto map = AMO::AdaptiveMapOptimizer(disk, target_disk);

    scale_problem(mesh, target, 1e0);

    auto map = AMO::AdaptiveMapOptimizer(mesh, target, AMO::QuadraticPenaltyTerm, AMO::Mine);

    std::cout << "running pipeline now" << std::endl;
    map.run_pipeline();



    //try with high resolution
    // TM lule_mesh;
    // OM::IO::read_mesh(lule_mesh, "frames/05reference_mesh.om");
    // auto lule_map = AMO::AdaptiveMapOptimizer(lule_mesh, target);

    // lule_map.one_step();



    //minimum_example();

    //classify_disk();
    //classify_rectangle();

    //collapse_test();

    //intuition_test();

    //split_test();





}





