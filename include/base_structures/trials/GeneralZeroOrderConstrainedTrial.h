#ifndef GENERAL_ZERO_ORDER_CONSTRAINED_TRIAL_H_
#define GENERAL_ZERO_ORDER_CONSTRAINED_TRIAL_H_

#include <vector>

#include <base_structures/trials/GeneralZeroOrderTrial.h>

namespace opt {
    struct GeneralZeroOrderConstrainedTrial : public GeneralZeroOrderTrial {
        std::vector<double> g;

        GeneralZeroOrderConstrainedTrial(double _x = 0.0, double _z = 0.0,
                                         const std::vector<double> &_g = std::vector<double>{})
            : GeneralZeroOrderTrial(_x, _z), g(_g) {};
    };
}

#endif // GENERAL_ZERO_ORDER_CONSTRAINED_TRIAL_H_
