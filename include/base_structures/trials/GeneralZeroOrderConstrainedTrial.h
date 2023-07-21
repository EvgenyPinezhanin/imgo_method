#ifndef GENERAL_ZERO_ORDER_CONSTRAINED_TRIAL_H_
#define GENERAL_ZERO_ORDER_CONSTRAINED_TRIAL_H_

#include <vector>

#include <base_structures/trials/GeneralZeroOrderTrial.h>

using std::vector;

namespace om {
    struct GeneralZeroOrderConstrainedTrial : public GeneralZeroOrderTrial {
        vector<double> g;

        GeneralZeroOrderConstrainedTrial(double _x = 0.0, double _z = 0.0, const vector<double> &_g)
                                        : GeneralZeroOrderTrial(_x, _z), g(_g) {};
    };
}

#endif // GENERAL_ZERO_ORDER_CONSTRAINED_TRIAL_H_
