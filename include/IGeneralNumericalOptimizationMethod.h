#ifndef I_GENERAL_NUMERICAL_OPTIMIZATION_METHOD_H_
#define I_GENERAL_NUMERICAL_OPTIMIZATION_METHOD_H_

#include <vector>

#include <IGeneralOptimizationMethod.h>

using std::vector;

template <typename Trial>
class IGeneralNumericalOptimizationMethod : public IGeneralOptimizationMethod {
protected:
    vector<Trial> trialPoints;

    double accuracy;
    int maxTrials;
    int maxFevals;

    virtual Trial newTrial(double x) = 0;
    virtual double newPoint() = 0;
    virtual double selectNewPoint() = 0;
    virtual bool stopConditions() = 0;

public:
    IGeneralNumericalOptimizationMethod(int _dimension, const vector<double> &_lowerBound, const vector<double> &_upBound,
                                 double _accuracy, int _maxTrials, int _maxFevals) : IGeneralOptimizationMethod(_dimension,
                                 _lowerBound, _upBound), trialPoints(0), accuracy(_accuracy), maxTrials(_maxTrials),
                                 maxFevals(_maxFevals) {};

    void setAccuracy(double _accuracy) { accuracy = _accuracy; };
    void setMaxTrials(int _maxTrials) { maxTrials = _maxTrials; };
    void setMaxFevals(int _maxFevals) { maxFevals = _maxFevals; };

    double getAccuracy() const { return accuracy; };
    int getMaxTrials() const { return maxTrials; };
    int getMaxFevals() const { return maxFevals; };

    void getTrialPoints(vector<Trial> &_trialPoints) const { _trialPoints = trialPoints; };
    int getNumberTrialPoints() const { return trialPoints.size(); };

    virtual void solve(int &numberTrials, int &numberFevals, vector<double> &X) = 0;
    virtual bool solveTest(int &numberTrials, int &numberFevals, vector<double> XOpt) = 0;
};

#endif // I_GENERAL_NUMERICAL_OPTIMIZATION_METHOD_H_
