#ifndef CONSTRAINED_TRIAL_H_
#define CONSTRAINED_TRIAL_H_

#include <Trials/Trial.h>

struct ConstrainedTrial : public Trial {
    int nu;

    ConstrainedTrial(double _x = 0.0, double _z = 0.0, int _nu = 0) : Trial(_x, _z), nu(_nu) {}
};

#endif // CONSTRAINED_TRIAL_H_
