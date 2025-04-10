// RemesherValentin.hh
#ifndef REMESHER_VALENTIN_HH
#define REMESHER_VALENTIN_HH

#include "RemesherBase.hh"

namespace AdaptiveMapOptimizer {

class RemesherValentin : public RemesherBase<RemesherValentin> {
public:
    using RemesherBase::RemesherBase;

    void split_needles_and_collapse_spikes(
        int& split_count,
        int& collapse_count,
        bool splits_allowed,
        bool collapses_allowed);


private:
    double split_angle_thresh = 130; // thresholds
    double collapse_side_len_ratio = 1. / 8.;

};


} // namespace AdaptiveMapOptimizer

#endif // REMESHER_VALENTIN_HH
