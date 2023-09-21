#ifndef PIYAVSKY_METHOD_H_
#define PIYAVSKY_METHOD_H_

#include <base_classes/opt_methods/IConstantEstimationOptMethod.h>
#include <opt_methods/ScanningMethod.h>

class PiyavskyMethod : public ScanningMethod, public opt::IConstantEstimationOptMethod {
public:
    struct Parameters : public ScanningMethod::Parameters {
        double reliability;

        Parameters(double _accuracy = 0.001, double _error = 0.001,
                   int _maxTrials = 1000, int _maxFevals = 1000,
                   double _reliability = 2.0):
            ScanningMethod::Parameters(_accuracy, _error, _maxTrials, _maxFevals),
            reliability(_reliability)
        {};
    };

protected:
    void estimatingConstant() override;

    void calcCharacteristic() override;

    double selectNewPoint() override;

public:
    PiyavskyMethod(const OneDimensionalProblem &_problem = OneDimensionalProblem(),
                   const Parameters &parameters = Parameters()):
        ScanningMethod(_problem, parameters),
        opt::IConstantEstimationOptMethod(parameters.reliability)
    {};

    void setParameters(const GeneralMethod::Parameters &parameters) override {
        ScanningMethod::setParameters(parameters);
        reliability = static_cast<const Parameters&>(parameters).reliability;
    };
    void getParameters(GeneralMethod::Parameters &parameters) const override {
        ScanningMethod::getParameters(parameters);
        static_cast<Parameters&>(parameters).reliability = reliability;
    };
};

#endif // PIYAVSKY_METHOD_H_
