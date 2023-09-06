#ifndef SCANNING_METHOD_H_
#define SCANNING_METHOD_H_

#include <base_classes/opt_methods/ICharacteristicOptMethod.h>
#include <base_classes/opt_methods/IGeneralNumericalOptMethod.h>
#include <base_structures/GeneralParametersNumericalOptMethod.h>
#include <result_methods/ResultMethod.h>
#include <opt_problems/OneDimensionalProblem.h>
#include <trials/Trial.h>

struct ScanningParameters : public opt::GeneralParametersNumericalOptMethod {
    ScanningParameters(double _accuracy = 0.001, double _error = 0.001,
                       int _maxTrials = 1000, int _maxFevals = 1000):
        opt::GeneralParametersNumericalOptMethod(_accuracy, _error, _maxTrials, _maxFevals)
    {};
};

class ScanningMethod : public opt::ICharacteristicOptMethod<Trial>,
    public opt::IGeneralNumericalOptMethod<Trial, double,
                                           OneDimensionalProblem, ResultMethod> {
private:
    void insertInSorted(const Trial &trial) override;

    void calcCharacteristic() override;

    Trial newTrial(const double &x) override;
    double selectNewPoint() override;

    double estimateSolution(double &x) const override;

    bool stopConditions() override;
    bool stopConditionsTest() override;

    void setDataInResultMethod(ResultMethod &result) const;

public:
    ScanningMethod(const OneDimensionalProblem &_problem = OneDimensionalProblem(),
                   const ScanningParameters &parameters = ScanningParameters()):
        opt::ICharacteristicOptMethod<Trial>(),
        opt::IGeneralNumericalOptMethod<Trial, double, OneDimensionalProblem,
                                        ResultMethod>(_problem, parameters)
    {};

    void setParameters(const opt::GeneralParametersNumericalOptMethod &parameters) override {
        opt::IGeneralNumericalOptMethod<Trial, double, OneDimensionalProblem,
                                        ResultMethod>::setParameters(parameters);
    };
    void getParameters(opt::GeneralParametersNumericalOptMethod &parameters) const override {
        opt::IGeneralNumericalOptMethod<Trial, double, OneDimensionalProblem,
                                        ResultMethod>::getParameters(parameters);
    };

    void solve(ResultMethod &result) override;
    bool solveTest(ResultMethod &result) override;
};

#endif // SCANNING_METHOD_H_
