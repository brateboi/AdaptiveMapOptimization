#ifndef HELPERS_HH
#define HELPERS_HH

#endif // HELPERS_HH

#include <includes.hh>

void minimum_example();
void energy_computation_test();


void collapse_test();
void constrain_non_original_vertices(TM &mesh_);
double total_area(TM &mesh_);
void intuition_test();
void split_test();



TM get_disk_mesh_without_interior_vertices(int n_vertices, double radius);
OptimizationTarget get_disk_optimization_target(TM &disk, double inner_radius, double outer_radius);

TM get_box_mesh_one_interior_vertex(int n_height_vertices, double height);
OptimizationTarget get_box_optimization_target(TM &box, double frequency, double amplitude);

void scale_problem(TM &mesh_, OptimizationTarget &target_, double scaling_factor);

std::vector<OM::HalfedgeHandle> is_collapse_okay(TM &mesh_, OM::EdgeHandle eh, double _epsilon = 1e-10);
bool triangle_flip_condition(TM &mesh_, OM::HalfedgeHandle &heh, double _epsilon);
