#ifndef STRESS_TESTS_HH
#define STRESS_TESTS_HH
#endif // STRESS_TESTS_HH

#include <includes.hh>

TM stress_test_disk(TM diskmesh_, double inner_radius, double outer_radius);
TM stress_test_box(TM boxmesh_, double frequency, double amplitude);


void do_stress_test_disk(double inner_radius, double outer_radius);
void do_stress_test_box(double frequency, double amplitude);

