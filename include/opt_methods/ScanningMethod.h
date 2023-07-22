#ifndef SCANNING_METHOD_H_
#define SCANNING_METHOD_H_

#include <base_classes/opt_methods/ICharacteristicOptMethod.h>
#include <base_classes/opt_methods/IGeneralNumericalOptMethod.h>
#include <result_methods/ResultMethod.h>
#include <opt_problems/OneDimensionalProblem.h>
#include <trials/Trial.h>

class ScanningMethod : public om::ICharacteristicOptMethod<Trial>,
    public om::IGeneralNumericalOptMethod<Trial, double, OneDimensionalProblem, ResultMethod> {
private:
    void insertInSorted(const Trial &trial) override;

    void calcCharacteristic() override;

    Trial newTrial(const double &x) override;
    double selectNewPoint() override;

    double estimateSolution(double &x) const override;

    bool stopConditions() override;
    bool stopConditionsTest() override;

    void setDataInResultMethod(ResultMethod &result);

public:
    ScanningMethod(const OneDimensionalProblem &_problem = OneDimensionalProblem(), double _accuracy = 0.001,
                   double _error = 0.001, int _maxTrials = 1000, int _maxFevals = 1000)
                  : om::ICharacteristicOptMethod<Trial>(),
                    om::IGeneralNumericalOptMethod<Trial, double, OneDimensionalProblem, ResultMethod>(_problem,
                    _accuracy, _error, _maxTrials, _maxFevals) {};

    void solve(ResultMethod &result) override;

    bool solveTest(ResultMethod &result) override;
};

#endif // SCANNING_METHOD_H_
