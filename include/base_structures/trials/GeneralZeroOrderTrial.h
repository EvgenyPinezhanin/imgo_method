#ifndef GENERAL_ZERO_ORDER_TRIAL_H_
#define GENERAL_ZERO_ORDER_TRIAL_H_

namespace om {
    struct GeneralZeroOrderTrial {
        double x, z;

        GeneralZeroOrderTrial(double _x = 0.0, double _z = 0.0) : x(_x), z(_z) {};
    };
}

#endif // GENERAL_ZERO_ORDER_TRIAL_H_
