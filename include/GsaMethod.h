#ifndef GSA_METHOD_H_
#define GSA_METHOD_H_

#include <vector>
#include <functional>

#include <IStronginOptimizationMethod.h>

using std::vector;
using std::function;

class GsaMethod : public IStronginOptimizationMethod<function<double(double)>> {
private:
    int insertInSorted(vector<Trial> &trials, Trial trial) override;
    double searchMin() const override;

    void calcCharacteristic() override;

    Trial newTrial(double x) override;
    double newPoint() override;
    double selectNewPoint() override;
    bool stopConditions() override;

    double compute(vector<double> X) const override;

public:
    GsaMethod(function<double(double)> _objFunction = nullptr, double _a = 0.0, double _b = 10.0, double _r = 2.0,
              double _accuracy = 0.001, int _maxTrials = 1000, int _maxFevals = 1000)
              : IStronginOptimizationMethod(_objFunction, 1, vector<double>{ _a }, vector<double>{ _b }, _r,
              _accuracy, _maxTrials, _maxFevals) {};
    
    void setLowerBound(double _a) { IGeneralOptimizationMethod::setLowerBound(vector<double>{ _a }); };
    void setUpBound(double _b) { IGeneralOptimizationMethod::setUpBound(vector<double>{ _b }); };
    void setBounds(double _a, double _b) { IGeneralOptimizationMethod::setBounds(vector<double>{ _a }, vector<double>{ _b }); };

    double getA() const { return lowerBound[0]; };
    double getB() const { return upBound[0]; };

    void solve(ResultMethod *result) override;

    bool solveTest(double xOpt, ResultMethod *result);
    bool solveTest(vector<double> XOpt, ResultMethod *result) override;
};

#endif // GSA_METHOD_H_
