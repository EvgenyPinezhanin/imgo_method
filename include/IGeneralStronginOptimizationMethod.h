#ifndef I_GENERAL_STRONGIN_OPTIMIZATION_METHOD_H_
#define I_GENERAL_STRONGIN_OPTIMIZATION_METHOD_H_

#include <vector>

#include <IGeneralCharacteristicOptimizationMethod.h>

using std::vector;

template <typename Trial>
class IGeneralStronginOptimizationMethod : public IGeneralCharacteristicOptimizationMethod<Trial> {
protected:
    vector<Trial> lastTrials;
    vector<int> lastTrialsPos;

    virtual int insertInSorted(vector<Trial> &trials, Trial trial) = 0;
    virtual double searchMin() const = 0;

public:
    IGeneralStronginOptimizationMethod(int _dimension, const vector<double> &_lowerBound, const vector<double> &_upBound,
                                       double _accuracy, int _maxTrials, int _maxFevals) : IGeneralCharacteristicOptimizationMethod<Trial>(
                                       _dimension, _lowerBound, _upBound, _accuracy, _maxTrials, _maxFevals), lastTrials(0),
                                       lastTrialsPos(0) {};
};

#endif // I_GENERAL_STRONGIN_OPTIMIZATION_METHOD_H_
