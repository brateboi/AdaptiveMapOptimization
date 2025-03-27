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

