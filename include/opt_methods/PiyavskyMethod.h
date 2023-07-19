#ifndef PIYAVSKY_METHOD_H_
#define PIYAVSKY_METHOD_H_

#include <base_classes/opt_methods/IConstantEstimationOptMethod.h>
#include <opt_methods/ScanningMethod.h>
#include <trials/Trial.h>

class PiyavskyMethod : public ScanningMethod, public IConstantEstimationOptMethod {
private:
    void calcCharacteristic() override;

    double selectNewPoint() override;

public:
    PiyavskyMethod(const OneDimensionalProblem &_problem = OneDimensionalProblem(), double _reliability = 1.0,
                   double _accuracy = 0.001, double _error = 0.001, int _maxTrials = 1000, int _maxFevals = 1000)
                  : ScanningMethod(_problem, _accuracy, _error, _maxTrials, _maxFevals),
                  IConstantEstimationOptMethod(_reliability) {};
};

#endif // PIYAVSKY_METHOD_H_
