#ifndef I_GENERAL_OPT_PROBLEM_H_
#define I_GENERAL_OPT_PROBLEM_H_

#include <vector>

namespace opt {
    template <typename ObjectType, typename SearchAreaType, typename PointType>
    class IGeneralOptProblem {
    protected:
        ObjectType object;
        SearchAreaType area;
        std::vector<PointType> optimalPoints;
        double optimalValue;

    public:
        using Point = PointType;

        IGeneralOptProblem(const ObjectType &_object, const SearchAreaType &_area,
                           const std::vector<PointType> &_optimalPoints, double _optimalValue):
            object(_object),
            area(_area),
            optimalPoints(_optimalPoints),
            optimalValue(_optimalValue)
        {};

        void setObject(const ObjectType &_object) { object = _object; };
        ObjectType getObject() const { return object; };

        SearchAreaType getSearchArea() const { return area; };
        void getOptimalPoints(std::vector<PointType> &_optimalPoints) const { _optimalPoints = optimalPoints; };
        double getOptimalValue() const { return optimalValue; };

        virtual double computeObjFunction(const PointType &x) const = 0;
    };
}

#endif // I_GENERAL_OPT_PROBLEM_H_
