#ifndef INDEX_TRIAL_H_
#define INDEX_TRIAL_H_

#include <base_structures/trials/GeneralZeroOrderTrial.h>

struct IndexTrial : public GeneralZeroOrderTrial {
    int nu;

    IndexTrial(double _x = 0.0, double _z = 0.0, int _nu = 0) : GeneralZeroOrderTrial(_x, _z), nu(_nu) {}
};

#endif // INDEX_TRIAL_H_
