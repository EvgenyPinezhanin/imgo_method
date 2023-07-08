#ifndef I_GENERAL_OPT_TASK_H_
#define I_GENERAL_OPT_TASK_H_

template <typename ObjectiveFunctionType, typename SearchAreaType, typename OptimalPointType>
class IGeneralOptTask {
protected:
    ObjectiveFunctionType objFunction;
    SearchAreaType area;
    OptimalPointType optimalPoint;

public:
    IGeneralOptTask(const ObjectiveFunctionType &_objFunction, const SearchAreaType &_area,
                    const OptimalPointType &_optPoint)
                    : objFunction(_objFunction), area(_area), optPoint(_optPoint) {};

    void setObjFunction(ObjectiveFunctionType _objFunction) { objFunction = _objFunction; };
    ObjectiveFunctionType getObjFunction() const { return objFunction; };

    void setSearchArea(const SearchAreaType &_area) { area = _area; };
    SearchAreaType getSearchArea() const { return area; };

    void setOptimalPoint(const OptimalPointType &_optPoint) { optPoint = _optPoint; };
    OptimalPointType getOptimalPoint() const { return optPoint; };

    virtual double computeObjFunction(double x) const = 0;
};

#endif // I_GENERAL_OPT_TASK_H_
