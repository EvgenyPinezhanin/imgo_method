#ifndef _INDEX_TRIAL_H_
#define _INDEX_TRIAL_H_

#include <general/structures/trials/GeneralZeroOrderTrial.h>

namespace opt {
    struct IndexTrial : public GeneralZeroOrderTrial {
        int nu;

        IndexTrial(double _x = 0.0, double _z = 0.0, int _nu = 0) : GeneralZeroOrderTrial(_x, _z), nu(_nu) {}
    };
}

#endif // _INDEX_TRIAL_H_
