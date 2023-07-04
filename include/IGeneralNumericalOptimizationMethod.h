#ifndef I_GENERAL_NUMERICAL_OPTIMIZATION_METHOD_H_
#define I_GENERAL_NUMERICAL_OPTIMIZATION_METHOD_H_

#include <vector>

#include <IGeneralOptimizationMethod.h>

using std::vector;

template <typename TaskOptimizationMethodType, typename TrialType>
class IGeneralNumericalOptimizationMethod : public IGeneralOptimizationMethod<TaskOptimizationMethodType> {
protected:
    vector<TrialType> trialPoints;

    double accuracy;
    int maxTrials, numberTrials;

    virtual TrialType newTrial(double x) = 0;
    virtual double newPoint() = 0;
    virtual double selectNewPoint() = 0;
    virtual bool stopConditions() = 0;

public:
    IGeneralNumericalOptimizationMethod(const TaskOptimizationMethodType &_task, double _accuracy, int _maxTrials, int _maxFevals)
                                        : IGeneralOptimizationMethod<TaskOptimizationMethodType>(_task, _maxFevals),
                                        trialPoints(0), accuracy(_accuracy), maxTrials(_maxTrials), numberTrials(0) {};

    void setAccuracy(double _accuracy) { accuracy = _accuracy; };
    double getAccuracy() const { return accuracy; };

    void setMaxTrials(int _maxTrials) { maxTrials = _maxTrials; };
    int getMaxTrials() const { return maxTrials; };

    void getTrialPoints(vector<TrialType> &_trialPoints) const { _trialPoints = trialPoints; };
    int getNumberTrialPoints() const { return trialPoints.size(); };
};

#endif // I_GENERAL_NUMERICAL_OPTIMIZATION_METHOD_H_
