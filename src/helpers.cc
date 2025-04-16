#include <adaptive_map.hh>


void constrain_non_original_vertices(TM &mesh_){
    OM::VPropHandleT<bool> fixed_prop;

    mesh_.add_property(fixed_prop, "fixed");

    mesh_.property(fixed_prop).set_persistent(true);

    // constrain the boundary as fixed in first iteration
    for (auto vh : mesh_.vertices()){
        mesh_.property(fixed_prop, vh) = mesh_.is_boundary(vh);
    }

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


std::vector<OM::HalfedgeHandle> is_collapse_okay(TM &mesh_, OM::EdgeHandle eh, double _epsilon){

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
