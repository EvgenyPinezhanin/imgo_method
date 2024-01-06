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

        IGeneralOptProblem(const ObjectType &_object, const SearchAreaType &_area = SearchAreaType(),
                           const std::vector<PointType> &_optimalPoints = std::vector<PointType>{},
                           double _optimalValue = 0.0)
            : object(_object), area(_area), optimalPoints(_optimalPoints), optimalValue(_optimalValue) {};

        void setObject(ObjectType &_object) { object = _object; };
        void getObject(ObjectType &_object) const { _object = object; };

        virtual SearchAreaType getSearchArea() const { return area; };
        virtual void getOptimalPoints(std::vector<PointType> &_optimalPoints) const { _optimalPoints = optimalPoints; };
        virtual double getOptimalValue() const { return optimalValue; };

        virtual double computeObjFunction(const PointType &x) const = 0;
    };
}

#endif // I_GENERAL_OPT_PROBLEM_H_
