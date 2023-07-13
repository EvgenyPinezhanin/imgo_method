#ifndef I_GENERAL_OPT_PROBLEM_H_
#define I_GENERAL_OPT_PROBLEM_H_

template <typename ObjectiveFunctionType, typename SearchAreaType, typename OptimalPointType>
class IGeneralOptProblem {
protected:
    ObjectiveFunctionType objFunction;
    SearchAreaType area;
    OptimalPointType optimalPoint;

public:
    IGeneralOptProblem(const ObjectiveFunctionType &_objFunction, const SearchAreaType &_area,
                       const OptimalPointType &_optimalPoint)
                      : objFunction(_objFunction), area(_area), optimalPoint(_optimalPoint) {};

    void setObjFunction(ObjectiveFunctionType _objFunction) { objFunction = _objFunction; };
    ObjectiveFunctionType getObjFunction() const { return objFunction; };

    void setSearchArea(const SearchAreaType &_area) { area = _area; };
    SearchAreaType getSearchArea() const { return area; };

    void setOptimalPoint(const OptimalPointType &_optimalPoint) { optimalPoint = _optimalPoint; };
    OptimalPointType getOptimalPoint() const { return optimalPoint; };

    virtual double computeObjFunction(double x) const = 0;
};

#endif // I_GENERAL_OPT_PROBLEM_H_
