#ifndef I_GENERAL_NUMERICAL_OPT_METHOD_H_
#define I_GENERAL_NUMERICAL_OPT_METHOD_H_

#include <vector>

#include <base_classes/opt_methods/IGeneralOptMethod.h>

using std::vector;

template <typename SolutionType, typename TrialType, typename OptTaskType, typename ResultMethodType, typename PointType>
class IGeneralNumericalOptMethod : public IGeneralOptMethod<OptTaskType, ResultMethodType> {
protected:
    vector<TrialType> trialPoints;

    double accuracy, error;
    int maxTrials, numberTrials;
    int maxFevals, numberFevals;

    virtual TrialType newTrial(const PointType &x) = 0;
    virtual PointType selectNewPoint() = 0;

    virtual SolutionType estimateSolution() const = 0;

    virtual bool stopConditions() = 0;
    virtual bool stopConditionsTest() = 0;

public:
    IGeneralNumericalOptMethod(const OptTaskType &_task, double _accuracy, double _error, int _maxTrials, int _maxFevals)
                              : IGeneralOptMethod<OptTaskType, ResultMethodType>(_task, _maxFevals), trialPoints(0),
                              accuracy(_accuracy), error(_error), maxTrials(_maxTrials), numberTrials(0),
                              maxFevals(_maxFevals), numberFevals(0) {};

    void setAccuracy(double _accuracy) { accuracy = _accuracy; };
    double getAccuracy() const { return accuracy; };

    void setError(double _error) { error = _error; };
    double getError() const { return error; };

    void setMaxTrials(int _maxTrials) { maxTrials = _maxTrials; };
    int getMaxTrials() const { return maxTrials; };

    void setMaxFevals(int _maxFevals) { maxFevals = _maxFevals; };
    int getMaxFevals() const { return maxFevals; };

    void getTrialPoints(vector<TrialType> &_trialPoints) const { _trialPoints = trialPoints; };
    int getNumberTrialPoints() const { return trialPoints.size(); };
};

#endif // I_GENERAL_NUMERICAL_OPT_METHOD_H_
