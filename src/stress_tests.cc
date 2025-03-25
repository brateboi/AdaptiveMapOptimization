#include <adaptive_map.hh>


TM stress_test_disk(TM diskmesh_, double inner_radius, double outer_radius) {

    using MyElem = SDEB2D_PH_AD;

    std::cout<<" ---------------------------------------------------------------------------------------- "<<std::endl;
    std::cout<<" *** trying to deform disk with vertices to inner radius " << inner_radius << " and outer radius " <<  outer_radius << std::endl;


    std::vector<Vec2d> P_ref(3);
    P_ref[0] << 0,0;
    P_ref[1] << 1,0;
    P_ref[2] << 0.5, 0.5 * std::sqrt(3.0);

    double scale_factor = 1;

    for(int i(0); i<3; i++){
        P_ref[i] *= scale_factor;
    }


    //low weight so the pull has the most influence and this only really serves as a barrier
    const double Dirichlet_w(1), Dirichlet_barrier_w(0), Dirichlet_barrier_min(1e8);
    const int max_inner_iters(100);
    const int n_vars(diskmesh_.n_vertices() * 2);

    COMISO::FiniteElementProblem fe_problem(n_vars);
    std::vector<COMISO::LinearConstraint> eq_constraints;


    COMISO::FiniteElementSet<MyElem> fe_barrier("barrier elements");

    //add the barrier elements
    for (auto fh: diskmesh_.faces()) {

        //std::cout<<" --------------------------------------------------- handling face "<<fh<<std::endl;

        double w = 1;
        double sizing_scale = 1.0; //compute_sizing_scale(ch);

        typename MyElem::VecI vi_barrier;
        typename MyElem::VecC vc_barrier;
        typename MyElem::VecV vx_barrier;

        //auto temp = mesh_.face(fh);
        std::vector<TM::VertexHandle> verts = std::vector(diskmesh_.fv_begin(fh), diskmesh_.fv_end(fh));

        //auto verts = mesh_.get_halfface_vertices(mesh_.halfface_handle(fh,0));


        for (int local_vid(0); local_vid < 3; local_vid++) {
            const auto &vh = verts[local_vid];

            int global_id = vh.idx() * 2;
            for (int i(0); i < 2; i++) {
                vi_barrier[2 * local_vid + i] = global_id + i;
            }

            for (int i(0); i < 2; i++) {
                fe_problem.x()[global_id + i] = diskmesh_.point(vh)[i];
                vx_barrier[2 * local_vid + i] = fe_problem.x()[global_id + i];
            }
        }

        auto setup_res = MyElem::compute_constants(Dirichlet_w, Dirichlet_barrier_w, Dirichlet_barrier_min,
                                                   P_ref,
                                                   vc_barrier);

        //std::cout<<" - constants vector: "<<vc_barrier.transpose()<<std::endl;

        if (setup_res) {
            std::cout << " ERROR: constants couldn't be computed for face: " <<fh<< std::endl;
            return TM();
        }

        fe_barrier.instances().add_element(vi_barrier,
                                           vc_barrier);
    }

    fe_problem.add_set(&fe_barrier);

    // find center
    OM::Vec3d center = OM::Vec3d(0,0,0);
    for (int i = 0; i<=22; ++i){
        auto vh = TM::VertexHandle(i);
        center += diskmesh_.point(vh);
    }
    center /= 23.0; // 23 vertices in total


    // constrain the fixed vertices, i.e. vertices 0-22 in alternating fashion to inner and outer radius from center.
    for (int i = 0; i<=22; ++i){
        auto vh = TM::VertexHandle(i);
        COMISO::LinearConstraint::SVectorNC coeffs0, coeffs1;
        coeffs0.resize(n_vars);
        coeffs1.resize(n_vars);

        int var_idx = vh.idx() * 2;
        coeffs0.coeffRef(var_idx + 0) = 1;
        coeffs1.coeffRef(var_idx + 1) = 1;

        auto pos = diskmesh_.point(vh);

        // pointing from center with length 1
        OM::Vec3d dir = pos - center;
        dir.normalize();

        if( i % 2 == 0) { // all even verts to inner radius
            OM::Vec3d new_pos = inner_radius * dir;

            eq_constraints.push_back(COMISO::LinearConstraint(coeffs0, -new_pos[0]));
            eq_constraints.push_back(COMISO::LinearConstraint(coeffs1, -new_pos[1]));
            std::cout << " - constrained inner vertex " << vh << " to " << new_pos << std::endl;
        } else { // all odd verts to outer radius
            OM::Vec3d new_pos = outer_radius * dir;

            eq_constraints.push_back(COMISO::LinearConstraint(coeffs0, -new_pos[0]));
            eq_constraints.push_back(COMISO::LinearConstraint(coeffs1, -new_pos[1]));
            std::cout << " - constrained outer vertex " << vh << " to " << new_pos << std::endl;
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

    for(auto v: diskmesh_.vertices()){
        int global_id = 2 * v.idx();

        typename OM::Vec3d new_pos;
        for(int i = 0; i<2; i++){
            new_pos[i] = fe_problem.x()[global_id + i];
        }
        diskmesh_.set_point(v, new_pos);
    }


    //return updated mesh positions
    return diskmesh_;


    std::cout<<" --------------------------------------------------------- "<<std::endl;

}

void do_stress_test_disk(double inner_radius, double outer_radius){
    TM diskmesh;
    OM::IO::read_mesh(diskmesh, "../../../meshes/disk.om" );

    TM deformed_disk = stress_test_disk(diskmesh, inner_radius, outer_radius);

    OM::IO::write_mesh(deformed_disk, "../../../meshes/disk_deformed.om");
}


TM stress_test_box(TM boxmesh_, double frequency, double amplitude) {

    using MyElem = SDEB2D_PH_AD;

    std::cout<<" ---------------------------------------------------------------------------------------- "<<std::endl;
    std::cout<<" *** trying to deform box with vertices to sine amplitude " << amplitude << std::endl;


    std::vector<Vec2d> P_ref(3);
    P_ref[0] << 0,0;
    P_ref[1] << 1,0;
    P_ref[2] << 0.5, 0.5 * std::sqrt(3.0);

    double scale_factor = 1;

    for(int i(0); i<3; i++){
        P_ref[i] *= scale_factor;
    }


    //low weight so the pull has the most influence and this only really serves as a barrier
    const double Dirichlet_w(1), Dirichlet_barrier_w(0), Dirichlet_barrier_min(1e1);
    const int max_inner_iters(100);
    const int n_vars(boxmesh_.n_vertices() * 2);

    COMISO::FiniteElementProblem fe_problem(n_vars);
    std::vector<COMISO::LinearConstraint> eq_constraints;


    COMISO::FiniteElementSet<MyElem> fe_barrier("barrier elements");

    //add the barrier elements
    for (auto fh: boxmesh_.faces()) {

        //std::cout<<" --------------------------------------------------- handling face "<<fh<<std::endl;

        double w = 1;
        double sizing_scale = 1.0; //compute_sizing_scale(ch);

        typename MyElem::VecI vi_barrier;
        typename MyElem::VecC vc_barrier;
        typename MyElem::VecV vx_barrier;

        //auto temp = mesh_.face(fh);
        std::vector<TM::VertexHandle> verts = std::vector(boxmesh_.fv_begin(fh), boxmesh_.fv_end(fh));

        //auto verts = mesh_.get_halfface_vertices(mesh_.halfface_handle(fh,0));


        for (int local_vid(0); local_vid < 3; local_vid++) {
            const auto &vh = verts[local_vid];

            int global_id = vh.idx() * 2;
            for (int i(0); i < 2; i++) {
                vi_barrier[2 * local_vid + i] = global_id + i;
            }

            for (int i(0); i < 2; i++) {
                fe_problem.x()[global_id + i] = boxmesh_.point(vh)[i];
                vx_barrier[2 * local_vid + i] = fe_problem.x()[global_id + i];
            }
        }

        auto setup_res = MyElem::compute_constants(Dirichlet_w, Dirichlet_barrier_w, Dirichlet_barrier_min,
                                                   P_ref,
                                                   vc_barrier);

        //std::cout<<" - constants vector: "<<vc_barrier.transpose()<<std::endl;

        if (setup_res) {
            std::cout << " ERROR: constants couldn't be computed for face: " <<fh<< std::endl;
            return TM();
        }

        fe_barrier.instances().add_element(vi_barrier,
                                           vc_barrier);
    }

    fe_problem.add_set(&fe_barrier);

    // sine deformation based on height
    auto deform_func = [&](double y){
        return amplitude * std::sin(frequency * y);
    };

    // constrain the fixed vertices, i.e. vertices 0-27 by using the y-cord as input to the sine-function as a deformation map.
    for (int i = 0; i<=27; ++i){
        auto vh = TM::VertexHandle(i);
        COMISO::LinearConstraint::SVectorNC coeffs0, coeffs1;
        coeffs0.resize(n_vars);
        coeffs1.resize(n_vars);

        int var_idx = vh.idx() * 2;
        coeffs0.coeffRef(var_idx + 0) = 1;
        coeffs1.coeffRef(var_idx + 1) = 1;

        auto pos = boxmesh_.point(vh);

        double x = pos[0];
        double y = pos[1];

        double new_x = x + deform_func(y);

        OM::Vec2d new_pos = OM::Vec2d(new_x, y);

        eq_constraints.push_back(COMISO::LinearConstraint(coeffs0, -new_pos[0]));
        eq_constraints.push_back(COMISO::LinearConstraint(coeffs1, -new_pos[1]));
        std::cout << " - constrained vertex " << vh << " to " << new_pos << std::endl;
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

    for(auto v: boxmesh_.vertices()){
        int global_id = 2 * v.idx();

        typename OM::Vec3d new_pos;
        for(int i = 0; i<2; i++){
            new_pos[i] = fe_problem.x()[global_id + i];
        }
        boxmesh_.set_point(v, new_pos);
    }


    //return updated mesh positions
    return boxmesh_;


    std::cout<<" --------------------------------------------------------- "<<std::endl;

}

void do_stress_test_box(double frequency, double amplitude){
    TM boxmesh;
    OM::IO::read_mesh(boxmesh, "../../../meshes/rectangle1.om" );

    TM deformed_box = stress_test_box(boxmesh, frequency, amplitude);

    OM::IO::write_mesh(deformed_box, "../../../meshes/rectangle1_deformed.om");
}
