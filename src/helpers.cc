#include <adaptive_map.hh>
#include <remesher.hh>




void constrain_non_original_vertices(TM &mesh_){
    OM::VPropHandleT<bool> fixed_prop;

    mesh_.add_property(fixed_prop, "fixed");

    mesh_.property(fixed_prop).set_persistent(true);

    // constrain the boundary as fixed in first iteration
    for (auto vh : mesh_.vertices()){
        mesh_.property(fixed_prop, vh) = mesh_.is_boundary(vh);
    }

}


TM optimize_loose_vertices(TM mesh_){
    using MyElem = SDEB2D_PH_AD;
    //using MyElem = COMISO::FiniteElementTinyAD<FoldoverFreeElement2D>;

    std::cout<<" --------------------------------------------------------- "<<std::endl;
    std::cout<<" *** optimizing loose vertices of mesh, i.e. all non-original or inner vertices"<<std::endl;


    std::vector<Vec2d> P_ref(3);
    P_ref[0] << 0,0;
    P_ref[1] << 1,0;
    P_ref[2] << 0.5, 0.5 * std::sqrt(3.0);

    // get average size of elements
    double average_size = total_area(mesh_) / mesh_.n_faces();
    double scale_factor = std::sqrt(average_size/0.433013);

    for(int i(0); i<3; i++){
        P_ref[i] *= scale_factor;
    }


    //low weight so the pull has the most influence and this only really serves as a barrier
    const double Dirichlet_w(1), Dirichlet_barrier_w(0), Dirichlet_barrier_min(1e8);
    const int max_inner_iters(100);
    const int n_vars(mesh_.n_vertices() * 2);

    COMISO::FiniteElementProblem fe_problem(n_vars);
    std::vector<COMISO::LinearConstraint> eq_constraints;


    COMISO::FiniteElementSet<MyElem> fe_barrier("barrier elements");

    //add the barrier elements
    for (auto fh: mesh_.faces()) {

        //std::cout<<" --------------------------------------------------- handling face "<<fh<<std::endl;

        double w = 1;
        double sizing_scale = 1.0; //compute_sizing_scale(ch);

        typename MyElem::VecI vi_barrier;
        typename MyElem::VecC vc_barrier;
        typename MyElem::VecV vx_barrier;

        //auto temp = mesh_.face(fh);
        std::vector<TM::VertexHandle> verts = std::vector(mesh_.fv_begin(fh), mesh_.fv_end(fh));

        //auto verts = mesh_.get_halfface_vertices(mesh_.halfface_handle(fh,0));


        for (int local_vid(0); local_vid < 3; local_vid++) {
            const auto &vh = verts[local_vid];

            int global_id = vh.idx() * 2;
            for (int i(0); i < 2; i++) {
                vi_barrier[2 * local_vid + i] = global_id + i;
            }

            for (int i(0); i < 2; i++) {
                fe_problem.x()[global_id + i] = mesh_.point(vh)[i];
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


    //constrain the fixed vertices
    for(auto v: mesh_.vertices()) {
        OM::VPropHandleT<bool> fixed_prop;
        if (mesh_.get_property_handle(fixed_prop, "fixed")){ // retrieving property fixed successful

            if (mesh_.property(fixed_prop, v)){
                COMISO::LinearConstraint::SVectorNC coeffs0, coeffs1;
                coeffs0.resize(n_vars);
                coeffs1.resize(n_vars);

                int var_idx = v.idx() * 2;
                coeffs0.coeffRef(var_idx + 0) = 1;
                coeffs1.coeffRef(var_idx + 1) = 1;

                auto pos = mesh_.point(v);


                eq_constraints.push_back(COMISO::LinearConstraint(coeffs0, -pos[0]));
                eq_constraints.push_back(COMISO::LinearConstraint(coeffs1, -pos[1]));
                std::cout << " - constrained vertex " << v << " to " << pos << std::endl;
            }

        } else {
            std::cerr << "Failed to retrieve property \" fixed \". Exiting... " << std::endl;
            exit(-1);
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

    for(auto v: mesh_.vertices()){
        int global_id = 2 * v.idx();

        typename OM::Vec3d new_pos;
        for(int i = 0; i<2; i++){
            new_pos[i] = fe_problem.x()[global_id + i];
        }
        mesh_.set_point(v, new_pos);
    }


    //return updated mesh positions
    return mesh_;


    std::cout<<" --------------------------------------------------------- "<<std::endl;

}


TM optimize_inner_vertices(TM mesh_) {

    using MyElem = SDEB2D_PH_AD;
    //using MyElem = COMISO::FiniteElementTinyAD<FoldoverFreeElement2D>;

    std::cout<<" --------------------------------------------------------- "<<std::endl;
    std::cout<<" *** optimizing inner vertices of mesh scaling the reference element by average area"<<std::endl;


    std::vector<Vec2d> P_ref(3);
    P_ref[0] << 0,0;
    P_ref[1] << 1,0;
    P_ref[2] << 0.5, 0.5 * std::sqrt(3.0);


    // get average size of elements
    //double average_size = total_area(mesh_) / mesh_.n_faces(); // 40 for the deformed rect case
    double average_size = 40.; // 40 for the deformed rect case
    double scale_factor = std::sqrt(average_size / 0.433013);

    std::cout << "SCALE FACTOR " << scale_factor << std::endl;

    for(int i(0); i<3; i++){
        P_ref[i] *= scale_factor;
    }


    //low weight so the pull has the most influence and this only really serves as a barrier
    const double Dirichlet_w(1), Dirichlet_barrier_w(0), Dirichlet_barrier_min(1e8);
    const int max_inner_iters(100);
    const int n_vars(mesh_.n_vertices() * 2);

    COMISO::FiniteElementProblem fe_problem(n_vars);
    std::vector<COMISO::LinearConstraint> eq_constraints;

    COMISO::FiniteElementSet<MyElem> fe_barrier("barrier elements");

    //add the barrier elements
    for (auto fh: mesh_.faces()) {

        //std::cout<<" --------------------------------------------------- handling face "<<fh<<std::endl;

        double w = 1;
        double sizing_scale = 1.0; //compute_sizing_scale(ch);

        typename MyElem::VecI vi_barrier;
        typename MyElem::VecC vc_barrier;
        typename MyElem::VecV vx_barrier;

        //auto temp = mesh_.face(fh);
        std::vector<TM::VertexHandle> verts = std::vector(mesh_.fv_begin(fh), mesh_.fv_end(fh));

        //auto verts = mesh_.get_halfface_vertices(mesh_.halfface_handle(fh,0));


        for (int local_vid(0); local_vid < 3; local_vid++) {
            const auto &vh = verts[local_vid];

            int global_id = vh.idx() * 2;
            for (int i(0); i < 2; i++) {
                vi_barrier[2 * local_vid + i] = global_id + i;
            }

            for (int i(0); i < 2; i++) {
                fe_problem.x()[global_id + i] = mesh_.point(vh)[i];
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
            return TM();
        }

        fe_barrier.instances().add_element(vi_barrier,
                                           vc_barrier);
    }

    fe_problem.add_set(&fe_barrier);


    //constrain the fixed vertices
    for(auto v: mesh_.vertices()) {
        if(mesh_.is_boundary(v)) {
            COMISO::LinearConstraint::SVectorNC coeffs0, coeffs1;
            coeffs0.resize(n_vars);
            coeffs1.resize(n_vars);

            int var_idx = v.idx() * 2;
            coeffs0.coeffRef(var_idx + 0) = 1;
            coeffs1.coeffRef(var_idx + 1) = 1;

            auto pos = mesh_.point(v);


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

    for(auto v: mesh_.vertices()){
        int global_id = 2 * v.idx();

        typename OM::Vec3d new_pos;
        for(int i = 0; i<2; i++){
            new_pos[i] = fe_problem.x()[global_id + i];
        }
        mesh_.set_point(v, new_pos);
    }


    //return updated mesh positions
    return mesh_;


    std::cout<<" --------------------------------------------------------- "<<std::endl;



}


void cleanup_disk_mesh(){

    TM mesh;
    OM::IO::read_mesh(mesh, "../../../meshes/disk.om" );
    auto heh = mesh.find_halfedge(OM::VertexHandle(21), OM::VertexHandle(22));


    OM::IO::write_mesh(mesh, "../../../meshes/cleaned_disk.om");
}

void classify_vertices(TM &mesh){
    OM::VPropHandleT<int> classification_prop;

    mesh.add_property(classification_prop, "classification");

    mesh.property(classification_prop).set_persistent(true);

    for (auto v : mesh.vertices()){
        if (mesh.is_boundary(v)){
            mesh.property(classification_prop, v) = 1;
        } else {
            mesh.property(classification_prop, v) = 0;
        }
    }

}

void classify_disk(){
    TM mesh;
    OM::IO::read_mesh(mesh, "../../../meshes/disk.om" );

    classify_vertices(mesh);

    OM::IO::write_mesh(mesh, "../../../meshes/disk_classified.om");
}

void classify_rectangle(){
    TM mesh;
    OM::IO::read_mesh(mesh, "../../../meshes/rectangle1.om" );

    classify_vertices(mesh);

    OM::IO::write_mesh(mesh, "../../../meshes/rectangle_classified.om");
}

TM get_minimum_mesh(){
    TM mesh;

    auto vh0 = mesh.add_vertex(OM::Vec3d(0,0,0));
    auto vh1 = mesh.add_vertex(OM::Vec3d(1,0,0));
    auto vh2 = mesh.add_vertex(OM::Vec3d(0.5,0.5*std::sqrt(3),0));
    auto vh3 = mesh.add_vertex(OM::Vec3d(0.5,0.4*std::sqrt(3),0));

    auto fh1 = mesh.add_face(vh0, vh1, vh3);
    auto fh2 = mesh.add_face(vh1, vh2, vh3);
    auto fh3 = mesh.add_face(vh2, vh0, vh3);

    std::cout << "Minimum mesh with #Vhs: " << mesh.n_vertices() << ", #Fhs: " << mesh.n_faces() << ", #Ehs: " << mesh.n_edges() << std::endl;

    for (auto v : mesh.vertices()){
        if (mesh.is_boundary(v)){
            std::cout << "found boundary vertex" << std::endl;
        }
    }

    return mesh;
}

double total_area(TM &mesh_){

    double area = 0;
    for (auto &f : mesh_.all_faces()){
        area += mesh_.calc_face_area(f);
    }

    return area;

}


double calculate_energy_of_mesh(TM &mesh_){

    using MyElem = SDEB2D_PH_AD;

    std::cout<<" --------------------------------------------------------- "<<std::endl;
    std::cout<<" *** Calculating energy with CoMISo"<<std::endl;


    std::vector<Vec2d> P_ref(3);
    P_ref[0] << 0,0;
    P_ref[1] << 1,0;
    P_ref[2] << 0.5, 0.5 * std::sqrt(3.0);

    // get average size of elements
    double average_size = total_area(mesh_) / mesh_.n_faces();
    std::cout << "area" << total_area(mesh_) << std::endl;

    // average_size is target_volume

    double scale_factor = std::sqrt(average_size / 0.433013);


    for(int i(0); i<3; i++){
        P_ref[i] *= scale_factor;
    }


    //low weight so the pull has the most influence and this only really serves as a barrier
    const double Dirichlet_w(1), Dirichlet_barrier_w(0), Dirichlet_barrier_min(1e8);
    const int max_inner_iters(100);
    const int n_vars(mesh_.n_vertices() * 2);

    COMISO::FiniteElementProblem fe_problem(n_vars);
    std::vector<COMISO::LinearConstraint> eq_constraints;


    COMISO::FiniteElementSet<MyElem> fe_barrier("barrier elements");

    //add the barrier elements
    for (auto fh: mesh_.faces()) {

        //std::cout<<" --------------------------------------------------- handling face "<<fh<<std::endl;

        double w = 1;
        double sizing_scale = 1.0; //compute_sizing_scale(ch);

        typename MyElem::VecI vi_barrier;
        typename MyElem::VecC vc_barrier;
        typename MyElem::VecV vx_barrier;

        //auto temp = mesh_.face(fh);
        std::vector<TM::VertexHandle> verts = std::vector(mesh_.fv_begin(fh), mesh_.fv_end(fh));

        //auto verts = mesh_.get_halfface_vertices(mesh_.halfface_handle(fh,0));


        for (int local_vid(0); local_vid < 3; local_vid++) {
            const auto &vh = verts[local_vid];

            int global_id = vh.idx() * 2;
            for (int i(0); i < 2; i++) {
                vi_barrier[2 * local_vid + i] = global_id + i;
            }

            for (int i(0); i < 2; i++) {
                fe_problem.x()[global_id + i] = mesh_.point(vh)[i];
                vx_barrier[2 * local_vid + i] = fe_problem.x()[global_id + i];
            }
        }

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

    std::cout<<" energies before optimization: "<<std::endl;
    fe_problem.print_objectives();


    std::cout<<" --------------------------------------------------------- "<<std::endl;


    auto SDE_Element = [&](OM::FaceHandle face) -> double {
        std::vector<OM::VertexHandle> verts;
        for (auto v_it = mesh_.fv_begin(face); v_it.is_valid(); ++v_it){
            verts.push_back(*v_it);
        }
        assert(verts.size() == 3);

        Vec2d v0 = Vec2d(mesh_.point(verts[0]).data());
        Vec2d v1 = Vec2d(mesh_.point(verts[1]).data());
        Vec2d v2 = Vec2d(mesh_.point(verts[2]).data());

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

    OM::FPropHandleT<double> energy_prop;
    mesh_.add_property(energy_prop, "energy");
    mesh_.property(energy_prop).set_persistent(true);


    double total_energy = 0;
    for (auto face : mesh_.faces()){
        double energy = SDE_Element(face);
        mesh_.property(energy_prop, face) = energy;


        total_energy += energy;
    }

    std::cout << "+++++ ENERGY CALCULATED MANUALLY: " << total_energy << std::endl;
    return total_energy;
}

void energy_computation_test(){
    TM mesh;

    auto vh0 = mesh.add_vertex(OM::Vec3d(0,0,0));
    auto vh1 = mesh.add_vertex(OM::Vec3d(2,0,0));
    auto vh2 = mesh.add_vertex(OM::Vec3d(0.5,0.1*std::sqrt(3),0));

    auto fh1 = mesh.add_face(vh0, vh1, vh2);

    double energy = calculate_energy_of_mesh(mesh);

    std::cout << "+++++ ENERGY CALCULATED MANUALLY: " << energy << std::endl;
}

void collapse_test(){
    TM mesh;
    mesh.request_edge_status();
    mesh.request_halfedge_status();
    mesh.request_vertex_status();
    mesh.request_face_status();
    constrain_non_original_vertices(mesh);

    auto a = mesh.add_vertex(OM::Vec3d(-2.2,-1,0));
    auto b = mesh.add_vertex(OM::Vec3d(-2,-1,0));
    auto c = mesh.add_vertex(OM::Vec3d(-2.48,-6.04,0));
    auto d = mesh.add_vertex(OM::Vec3d(-5.04,-2.82,0));
    auto e = mesh.add_vertex(OM::Vec3d(-2.5,0.7,0));
    auto f = mesh.add_vertex(OM::Vec3d(-0.4,0.78,0));
    auto g = mesh.add_vertex(OM::Vec3d(0.18,-2.84,0));

    auto fh1 = mesh.add_face(a, c, b);
    auto fh2 = mesh.add_face(a, d, c);
    auto fh3 = mesh.add_face(a, e, d);
    auto fh4 = mesh.add_face(a, b, e);
    auto fh5 = mesh.add_face(b, f, e);
    auto fh6 = mesh.add_face(b, g, f);
    auto fh7 = mesh.add_face(b, c, g);

    OM::IO::write_mesh(mesh, "../../../meshes/before_collapse_test.om");

    TM after_collapse = one_refine_step(mesh);

    after_collapse.garbage_collection();


    OM::IO::write_mesh(after_collapse, "../../../meshes/after_collapse_test.om");
}


void intuition_test(){
    TM before;
    TM after;

    OM::IO::read_mesh(before, "../../../meshes/before_collapse_test.om");
    OM::IO::read_mesh(after, "../../../meshes/after_collapse_test.om");

    calculate_energy_of_mesh(before);
    calculate_energy_of_mesh(after);

    TM before_optimized = optimize_inner_vertices(before);
    TM after_optimized = optimize_inner_vertices(after);

    calculate_energy_of_mesh(before_optimized);
    calculate_energy_of_mesh(after_optimized);

    OM::IO::write_mesh(after_optimized, "../../../meshes/optimized_after_collapse_test.om");
    OM::IO::write_mesh(before_optimized, "../../../meshes/optimized_before_collapse_test.om");

    OM::IO::write_mesh(after, "../../../meshes/after_collapse_test.om");
    OM::IO::write_mesh(before, "../../../meshes/before_collapse_test.om");

}

void split_test(){
    TM mesh;
    mesh.request_edge_status();
    mesh.request_halfedge_status();
    mesh.request_vertex_status();
    mesh.request_face_status();
    constrain_non_original_vertices(mesh);

    auto a = mesh.add_vertex(OM::Vec3d(-2.64,-0.18,0));
    auto b = mesh.add_vertex(OM::Vec3d(-0.1,0.18,0));
    auto c = mesh.add_vertex(OM::Vec3d(3.48,-0.32,0));
    auto d = mesh.add_vertex(OM::Vec3d(0.58,-6.16,0));
    auto e = mesh.add_vertex(OM::Vec3d(-8.46,-1.2,0));
    auto f = mesh.add_vertex(OM::Vec3d(-2.84,3.72,0));
    auto g = mesh.add_vertex(OM::Vec3d(2.92,3.68,0));
    auto h = mesh.add_vertex(OM::Vec3d(10,-1,0));

    auto fh1 = mesh.add_face(a, c, b);
    auto fh2 = mesh.add_face(a, b, f);
    auto fh3 = mesh.add_face(a, f, e);
    auto fh4 = mesh.add_face(a, e, d);
    auto fh5 = mesh.add_face(a, d, c);
    //auto fh6 = mesh.add_face(b, g, f); // not in one ring of edge
    auto fh7 = mesh.add_face(b, c, g);
    auto fh8 = mesh.add_face(c, d, h);
    auto fh9 = mesh.add_face(c, h, g);

    OM::IO::write_mesh(mesh, "../../../meshes/before_split_test.om");

    TM after_split = one_refine_step(mesh);


    after_split.garbage_collection();

    TM before_optimized = optimize_inner_vertices(mesh);
    TM after_optimized = optimize_inner_vertices(after_split);


    calculate_energy_of_mesh(mesh);
    calculate_energy_of_mesh(after_split);
    calculate_energy_of_mesh(before_optimized);
    calculate_energy_of_mesh(after_optimized);

    OM::IO::write_mesh(after_optimized, "../../../meshes/optimized_after_split_test.om");
    OM::IO::write_mesh(before_optimized, "../../../meshes/optimized_before_split_test.om");


    OM::IO::write_mesh(after_split, "../../../meshes/after_split_test.om");
}


void minimum_example(){

    TM mesh = get_minimum_mesh();

    OM::IO::write_mesh(mesh, "../../../meshes/base_mesh.om");

    TM optimized_mesh = optimize_inner_vertices(mesh);

    OM::IO::write_mesh(optimized_mesh, "../../../meshes/optimized_mesh.om");
}



TM get_disk_mesh_without_interior_vertices(int n_vertices, double radius){

    bool balanced_mode = true;

    TM disk_mesh = TM();

    auto get_pos = [&](int i) -> OM::Vec3d{
        double theta = ((double)i / (double)n_vertices) * 2* M_PI;
        return OM::Vec3d(radius * std::cos(theta), radius * std::sin(theta), 0);
    };

    int n_faces = n_vertices - 2;

    if (n_faces <= 0 ){
        std::cerr << "Need atleast 3 vertices" << std::endl;
        exit(-1);
    }

    std::vector<OM::VertexHandle> verts;
    for (int i = 0; i < n_vertices; i++){
        verts.push_back(disk_mesh.add_vertex(get_pos(i)));
    }

    if(balanced_mode){
        int passes = n_vertices / 2;
        if (n_vertices % 2  == 0){ // even --> top and bottom face, + left + right

            // top
            disk_mesh.add_face(verts[0], verts[1], verts[n_vertices-1]);

            // bottom
            disk_mesh.add_face(verts[passes], verts[passes+1], verts[passes-1]);

            for (int i = 1; i < passes-1; i++){

                disk_mesh.add_face(verts[i], verts[i+1], verts[n_vertices-i]);
                disk_mesh.add_face(verts[i+1], verts[n_vertices-i-1], verts[n_vertices-i]);
            }



        } else { // odd --> top, left + right
            // top
            disk_mesh.add_face(verts[0], verts[1], verts[n_vertices-1]);

            for (int i = 1; i < passes - 1; i++){
                disk_mesh.add_face(verts[i], verts[i+1], verts[n_vertices-i]);
                disk_mesh.add_face(verts[i+1], verts[n_vertices-i-1], verts[n_vertices-i]);
            }
        }
    } else {
        for (int i = 1; i < n_vertices - 1; i++){
            disk_mesh.add_face(verts[0], verts[i], verts[i+1]);
        }
    }

    return disk_mesh;
}

TM get_box_mesh_one_interior_vertex(int n_height_vertices, double height){
    TM box_mesh;

    double width = 1;

    OM::Vec3d center = OM::Vec3d(width/2, height/2, 0);

    for (int i = 0; i <= n_height_vertices; i++){

        OM::Vec3d pos1 = OM::Vec3d(0, (double)i/n_height_vertices * height, 0);
        OM::Vec3d pos2 = OM::Vec3d(1, (double)i/n_height_vertices * height, 0);

        box_mesh.add_vertex(pos1);
        box_mesh.add_vertex(pos2);
    }

    auto c_vh = box_mesh.add_vertex(center);

    for (int i = 0; i < 2*n_height_vertices; i+= 2){
        std::cout << " " << i << " " << c_vh.idx() << " " << i+2 << std::endl;
        std::cout << " " << i+1 << " " << i+3 << " " << c_vh.idx() <<  std::endl;


        box_mesh.add_face(OM::VertexHandle(i), OM::VertexHandle(c_vh), OM::VertexHandle(i+2));
        box_mesh.add_face(OM::VertexHandle(i+1), OM::VertexHandle(i+3), OM::VertexHandle(c_vh));
    }

    box_mesh.add_face(OM::VertexHandle(0), OM::VertexHandle(1), OM::VertexHandle(c_vh));
    box_mesh.add_face(OM::VertexHandle(2*n_height_vertices), OM::VertexHandle(c_vh), OM::VertexHandle(2*n_height_vertices + 1));


    return box_mesh;
}


OptimizationTarget get_disk_optimization_target(TM &disk, double inner_radius, double outer_radius){

    OptimizationTarget target;

    OM::Vec3d center = OM::Vec3d(0,0,0);

    for (int i = 0; i < disk.n_vertices(); i++){
        auto vh = OM::VertexHandle(i);
        auto pos = disk.point(vh);

        // pointing from center with length 1
        OM::Vec3d dir = pos - center;
        dir.normalize();

        OM::Vec3d new_pos;

        if( i % 2 == 0) { // all even verts to inner radius
            new_pos = inner_radius * dir;
        } else { // all odd verts to outer radius
            new_pos = outer_radius * dir;
        }
        target.push_back(std::pair(vh, new_pos));
    }
    return target;
}

OptimizationTarget get_box_optimization_target(TM &box, double frequency, double amplitude){
    // sine deformation based on height
    auto deform_func = [&](double y){
        return amplitude * std::sin(frequency * y);
    };

    OptimizationTarget target;


    // constrain the fixed vertices, i.e. vertices 0-27 by using the y-cord as input to the sine-function as a deformation map.
    for (int i = 0; i < box.n_vertices(); ++i){
        auto vh = TM::VertexHandle(i);


        auto pos = box.point(vh);

        double x = pos[0];
        double y = pos[1];

        double new_x = x + deform_func(y);

        OM::Vec3d new_pos = OM::Vec3d(new_x, y, 0);

        target.push_back(std::pair(vh, new_pos));
    }

    return target;
}


void scale_problem(TM &mesh_, OptimizationTarget &target_, double scaling_factor){
    for (const auto vh : mesh_.vertices()){
        OM::Vec3d new_pos = mesh_.point(vh) * scaling_factor;
        mesh_.set_point(vh, new_pos);
    }

    for (size_t i = 0; i < target_.size(); i++){
        auto pair = target_.at(i);
        OM::Vec3d new_pos = pair.second * scaling_factor;

        target_.at(i) = {pair.first, new_pos};

    }
}
