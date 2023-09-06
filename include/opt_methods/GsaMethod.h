#ifndef GSA_METHOD_H_
#define GSA_METHOD_H_

#include <opt_methods/PiyavskyMethod.h>

struct GsaParameters : public PiyavskyParameters {
    GsaParameters(double _accuracy = 0.001, double _error = 0.001,
                  int _maxTrials = 1000, int _maxFevals = 1000,
                  double _reliability = 2.0):
        PiyavskyParameters(_accuracy, _error, _maxTrials, _maxFevals, _reliability)
    {};
};


class GsaMethod : public PiyavskyMethod {
private:
    void calcCharacteristic() override;

public:
    GsaMethod(const OneDimensionalProblem &_problem = OneDimensionalProblem(),
              const GsaParameters &parameters = GsaParameters()):
        PiyavskyMethod(_problem, parameters)
    {};
};

#endif // GSA_METHOD_H_
