#ifndef REMESHERNAIVE_HH
#define REMESHERNAIVE_HH

#include "RemesherBase.hh"

namespace AdaptiveMapOptimizer {

class RemesherNaive : public RemesherBase<RemesherNaive> {
public:
    using RemesherBase::RemesherBase;

    void split_needles_and_collapse_spikes(
        int& split_count,
        int& collapse_count,
        bool splits_allowed,
        bool collapses_allowed);


private:

    double split_angle_thresh = 50;

};


} // namespace AdaptiveMapOptimizer

#endif // REMESHER_NAIVE_HH
