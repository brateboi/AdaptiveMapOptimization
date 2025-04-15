#ifndef HELPERS_HH
#define HELPERS_HH

#endif // HELPERS_HH

#include <includes.hh>

void minimum_example();
void energy_computation_test();
double calculate_energy_of_mesh(TM &mesh_);
TM optimize_inner_vertices(TM mesh_);
TM optimize_loose_vertices(TM mesh_);
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
