#ifndef ONE_DIMENSIONAL_OPTIMIZATION_TASK_H_
#define ONE_DIMENSIONAL_OPTIMIZATION_TASK_H_

template <typename ObjectiveFunctionType>
class OneDimensionalOptimizationTask {
protected:
    ObjectiveFunctionType objFunction;
    double lowerBound, upBound; // area of search

public:
    OneDimensionalOptimizationTask(ObjectiveFunctionType _objFunction, double _lowerBound, double _upBound)
                                   : objFunction(_objFunction), lowerBound(_lowerBound), upBound(_upBound) {};

    void setObjFunction(ObjectiveFunctionType _objFunction) { objFunction = _objFunction; };
    ObjectiveFunctionType getObjFunction() const { return objFunction; };

    void setLowerBound(double _lowerBound) { lowerBound = _lowerBound; };
    double getLowerBound() const { return lowerBound; };

    void setUpBound(double _upBound) { upBound = _upBound; };
    double getUpBound(double _upBound) const { return upBound; };

    void setBounds(double _lowerBound, double &_upBound) { lowerBound = _lowerBound; upBound = _upBound; };
    void getBounds(double &_lowerBound, double &_upBound) const { _lowerBound = lowerBound; _upBound = upBound; };

    virtual double computeObjFunction(double x) const = 0;
};

#endif // ONE_DIMENSIONAL_OPTIMIZATION_TASK_H_
