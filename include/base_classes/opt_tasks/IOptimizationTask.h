#ifndef I_OPTIMIZATION_TASK_H_
#define I_OPTIMIZATION_TASK_H_

template <typename ObjectiveFunctionType, typename SearchAreaType, typename OptimalPointType>
class IOptimizationTask {
protected:
    ObjectiveFunctionType objFunction;
    SearchAreaType area;
    OptimalPointType optPoint;

public:
    IOptimizationTask(const ObjectiveFunctionType &_objFunction, const SearchAreaType &_area, const OptimalPointType &_optPoint)
                      : objFunction(_objFunction), area(_area), optPoint(_optPoint) {};

    void setObjFunction(ObjectiveFunctionType _objFunction) { objFunction = _objFunction; };
    ObjectiveFunctionType getObjFunction() const { return objFunction; };

    void setSearchArea(const SearchAreaType &_area) { area = _area; };
    SearchAreaType getSearchArea() const { return area; };

    virtual double computeObjFunction(double x) const = 0;
};

#endif // I_OPTIMIZATION_TASK_H_
