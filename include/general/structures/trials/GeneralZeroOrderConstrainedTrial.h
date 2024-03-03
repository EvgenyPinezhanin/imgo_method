#ifndef _GENERAL_ZERO_ORDER_CONSTRAINED_TRIAL_H_
#define _GENERAL_ZERO_ORDER_CONSTRAINED_TRIAL_H_

#include <vector>

#include <general/structures/trials/GeneralZeroOrderTrial.h>

namespace opt {
    struct GeneralZeroOrderConstrainedTrial : public GeneralZeroOrderTrial {
        std::vector<double> g;

        GeneralZeroOrderConstrainedTrial(double _x = 0.0, double _z = 0.0,
                                         const std::vector<double> &_g = std::vector<double>{})
            : GeneralZeroOrderTrial(_x, _z), g(_g) {};
    };
}

#endif // _GENERAL_ZERO_ORDER_CONSTRAINED_TRIAL_H_
