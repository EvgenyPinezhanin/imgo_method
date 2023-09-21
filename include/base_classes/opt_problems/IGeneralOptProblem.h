#ifndef I_GENERAL_OPT_PROBLEM_H_
#define I_GENERAL_OPT_PROBLEM_H_

#include <vector>

using std::vector;

namespace opt {
    template <typename ObjectiveFunctionType, typename SearchAreaType, typename PointType>
    class IGeneralOptProblem {
    protected:
        ObjectiveFunctionType objFunction;
        SearchAreaType area;
        vector<PointType> optimalPoints;
        double optimalValue;

    public:
        using Point = PointType;

        IGeneralOptProblem(const ObjectiveFunctionType &_objFunction, const SearchAreaType &_area,
                           const vector<PointType> &_optimalPoints, double _optimalValue)
                          : objFunction(_objFunction), area(_area), optimalPoints(_optimalPoints),
                          optimalValue(_optimalValue) {};

        void setObjFunction(const ObjectiveFunctionType &_objFunction) { objFunction = _objFunction; };
        ObjectiveFunctionType getObjFunction() const { return objFunction; };

        void setSearchArea(const SearchAreaType &_area) { area = _area; };
        SearchAreaType getSearchArea() const { return area; };

        void setOptimalPoints(const vector<PointType> &_optimalPoints) { optimalPoints = _optimalPoints; };
        void getOptimalPoints(vector<PointType> &_optimalPoints) const { _optimalPoints = optimalPoints; };

        void setOptimalValue(double _optimalValue) { optimalValue = _optimalValue; };
        double getOptimalValue() const { return optimalValue; };

        virtual double computeObjFunction(const PointType &x) const = 0;
    };
}

#endif // I_GENERAL_OPT_PROBLEM_H_
