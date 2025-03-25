#include <includes.hh>
#include <adaptive_map.hh>
#include <stress_tests.hh>
#include <remesher.hh>

void hello_world_adopt()
{
    std::cout << "HELLO FROM MY LIBRARY" << std::endl;
}




int main(int argc, char const *argv[])
{
    std::cout << "Hello World " << std::endl;

    //energy_computation_test();

    //minimum_example();

    //classify_disk();
    //classify_rectangle();

    //cleanup_disk_mesh();

    //do_stress_test_disk(0.10, 5.0);

    // do_stress_test_box(2.0, 100.0);

    //collapse_test();

    //intuition_test();


    TM deformed_mesh;

    OM::IO::read_mesh(deformed_mesh, "../../../meshes/rectangle1_deformed.om" );

    calculate_energy_of_mesh(deformed_mesh);

    deformed_mesh.request_edge_status();
    deformed_mesh.request_halfedge_status();
    deformed_mesh.request_vertex_status();
    deformed_mesh.request_face_status();


    constrain_non_original_vertices(deformed_mesh);
    refine_and_optimize(deformed_mesh);


}





