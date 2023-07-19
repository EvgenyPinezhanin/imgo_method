#ifndef I_GENERAL_OPT_PROBLEM_H_
#define I_GENERAL_OPT_PROBLEM_H_

#include <vector>

using std::vector;

template <typename ObjectiveFunctionType, typename SearchAreaType, typename OptimalPointType>
class IGeneralOptProblem {
protected:
    ObjectiveFunctionType objFunction;
    SearchAreaType area;
    vector<OptimalPointType> optimalPoints;

public:
    IGeneralOptProblem(const ObjectiveFunctionType &_objFunction, const SearchAreaType &_area,
                       const vector<OptimalPointType> &_optimalPoints)
                      : objFunction(_objFunction), area(_area), optimalPoints(_optimalPoints) {};

    void setObjFunction(ObjectiveFunctionType _objFunction) { objFunction = _objFunction; };
    ObjectiveFunctionType getObjFunction() const { return objFunction; };

    void setSearchArea(const SearchAreaType &_area) { area = _area; };
    SearchAreaType getSearchArea() const { return area; };

    void setOptimalPoints(const vector<OptimalPointType> &_optimalPoints) { optimalPoints = _optimalPoints; };
    void getOptimalPoints(vector<OptimalPointType> &_optimalPoints) const { _optimalPoints = optimalPoints; };

    virtual double computeObjFunction(double x) const = 0;
};

#endif // I_GENERAL_OPT_PROBLEM_H_
