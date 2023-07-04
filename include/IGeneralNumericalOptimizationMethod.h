#ifndef I_GENERAL_NUMERICAL_OPTIMIZATION_METHOD_H_
#define I_GENERAL_NUMERICAL_OPTIMIZATION_METHOD_H_

#include <vector>

#include <IGeneralOptimizationMethod.h>

using std::vector;

template <typename ObjectiveFunction, typename Trial>
class IGeneralNumericalOptimizationMethod : public IGeneralOptimizationMethod<ObjectiveFunction> {
protected:
    vector<Trial> trialPoints;

    double accuracy;
    int maxTrials, numberTrials;

    virtual Trial newTrial(double x) = 0;
    virtual double newPoint() = 0;
    virtual double selectNewPoint() = 0;
    virtual bool stopConditions() = 0;

public:
    IGeneralNumericalOptimizationMethod(ObjectiveFunction _objFunction, int _dimension, const vector<double> &_lowerBound,
                                        const vector<double> &_upBound, double _accuracy, int _maxTrials, int _maxFevals)
                                        : IGeneralOptimizationMethod<ObjectiveFunction>(_objFunction, _dimension, _lowerBound,
                                        _upBound, _maxFevals), trialPoints(0), accuracy(_accuracy), maxTrials(_maxTrials),
                                        numberTrials(0) {};

    void setAccuracy(double _accuracy) { accuracy = _accuracy; };
    double getAccuracy() const { return accuracy; };

    void setMaxTrials(int _maxTrials) { maxTrials = _maxTrials; };
    int getMaxTrials() const { return maxTrials; };

    void getTrialPoints(vector<Trial> &_trialPoints) const { _trialPoints = trialPoints; };
    int getNumberTrialPoints() const { return trialPoints.size(); };
};

#endif // I_GENERAL_NUMERICAL_OPTIMIZATION_METHOD_H_
