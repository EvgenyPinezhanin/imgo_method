#ifndef GSA_METHOD_H_
#define GSA_METHOD_H_

#include <opt_methods/PiyavskyMethod.h>

class GsaMethod : public PiyavskyMethod {
private:
    void calcCharacteristic() override;

public:
    GsaMethod(const OneDimensionalProblem &_problem = OneDimensionalProblem(), double _reliability = 1.0,
              double _accuracy = 0.001, double _error = 0.001, int _maxTrials = 1000, int _maxFevals = 1000)
             : PiyavskyMethod(_problem, _reliability, _accuracy, _error, _maxTrials, _maxFevals) {};
};

#endif // GSA_METHOD_H_
