#ifndef PIYAVSKY_METHOD_H_
#define PIYAVSKY_METHOD_H_

#include <base_classes/opt_methods/IConstantEstimationOptMethod.h>
#include <opt_methods/ScanningMethod.h>
#include <trials/Trial.h>

struct PiyavskyParameters : public ScanningParameters {
    double reliability;

    PiyavskyParameters(double _accuracy = 0.001, double _error = 0.001,
                       int _maxTrials = 1000, int _maxFevals = 1000,
                       double _reliability = 2.0):
        ScanningParameters(_accuracy, _error, _maxTrials, _maxFevals),
        reliability(_reliability)
    {};
};

class PiyavskyMethod : public ScanningMethod, public opt::IConstantEstimationOptMethod {
private:
    void calcCharacteristic() override;

    double selectNewPoint() override;

public:
    PiyavskyMethod(const OneDimensionalProblem &_problem = OneDimensionalProblem(),
                   const PiyavskyParameters &parameters = PiyavskyParameters()):
        ScanningMethod(_problem, static_cast<ScanningParameters>(parameters)),
        opt::IConstantEstimationOptMethod(parameters.reliability)
    {};

    void setParameters(const opt::GeneralParametersNumericalOptMethod &parameters) override {
        ScanningMethod::setParameters(parameters);
        reliability = static_cast<const PiyavskyParameters&>(parameters).reliability;
    };
    void getParameters(opt::GeneralParametersNumericalOptMethod &parameters) const override {
        ScanningMethod::getParameters(parameters);
        static_cast<PiyavskyParameters&>(parameters).reliability = reliability;
    };
};

#endif // PIYAVSKY_METHOD_H_
