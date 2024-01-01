#ifndef INDEX_TRIAL_H_
#define INDEX_TRIAL_H_

#include <base_structures/trials/GeneralZeroOrderTrial.h>

namespace opt {
    struct IndexTrial : public GeneralZeroOrderTrial {
        size_t nu;

        IndexTrial(double _x = 0.0, double _z = 0.0, size_t _nu = 0) : GeneralZeroOrderTrial(_x, _z), nu(_nu) {}
    };
}

#endif // INDEX_TRIAL_H_
