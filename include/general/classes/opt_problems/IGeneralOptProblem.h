#ifndef _I_GENERAL_OPT_PROBLEM_H_
#define _I_GENERAL_OPT_PROBLEM_H_

#include <vector>

namespace opt {
    template <typename SearchAreaType, typename PointType>
    class IGeneralOptProblem {
    protected:
        SearchAreaType area;

    public:
        using Point = PointType;

        IGeneralOptProblem(const SearchAreaType &_area = SearchAreaType()) : area(_area) {};

        virtual SearchAreaType getSearchArea() const { return area; };
        virtual void getOptimalPoints(std::vector<PointType> &_optimalPoints) const = 0;
        virtual double getOptimalValue() const = 0;

        virtual double computeObjectiveFunction(const PointType &x) const = 0;
    };
}

#endif // _I_GENERAL_OPT_PROBLEM_H_
