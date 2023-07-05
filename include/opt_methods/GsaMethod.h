#ifndef GSA_METHOD_H_
#define GSA_METHOD_H_

#include <base_classes/opt_methods/IStronginOptimizationMethod.h>
#include <opt_tasks/OneDimensionalTask.h>
#include <ResultMethods/StronginResultMethod.h>
#include <Trials/Trial.h>

class GsaMethod : public IStronginOptimizationMethod<double, OneDimensionalTask, StronginResultMethod, double> {
private:
    int insertInSorted(Trial trial) override;
    double searchMin() const override;

    void calcCharacteristic() override;

    Trial newTrial(const double &x) override;
    double newPoint() override;
    double selectNewPoint() override;
    bool stopConditions() override;

public:
    GsaMethod(const OneDimensionalTask &_task, double _r = 2.0, double _accuracy = 0.001, int _maxTrials = 1000,
              int _maxFevals = 1000)
              : IStronginOptimizationMethod(_task, _r, _accuracy, _maxTrials, _maxFevals) {};

    void solve(StronginResultMethod &result) override;

    bool solveTest(const double &xOpt, StronginResultMethod &result) override;
};

#endif // GSA_METHOD_H_
