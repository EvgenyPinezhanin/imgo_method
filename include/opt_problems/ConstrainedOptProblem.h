#ifndef _CONSTRAINED_OPT_PROBLEM_H_
#define _CONSTRAINED_OPT_PROBLEM_H_

#include <vector>
#include <functional>

#include <opt_problems/base/IBaseConstrainedOptProblem.h>
#include <general/classes/opt_problems/IObjOptProblem.h>
#include <general/structures/search_areas/OneDimensionalSearchArea.h>
#include <general/structures/search_areas/MultiDimensionalSearchArea.h>

template <typename SearchAreaType, typename PointType>
class ConstrainedOptProblem:
    public IBaseConstrainedOptProblem<SearchAreaType, PointType>,
    public opt::IObjOptProblem<std::function<double(PointType, size_t)>>
{
public:
    ConstrainedOptProblem(const std::function<double(PointType, size_t)> &_object = nullptr,
                          const std::string &_blockName = "", size_t _problemNumber = 0, size_t _numberConstraints = 0,
                          const SearchAreaType &_area = SearchAreaType(),
                          const std::vector<PointType> &_optimalPoints = std::vector<PointType>{},
                          double _optimalValue = 0.0, double _objectiveLipschitzConstant = 0.0,
                          const std::vector<PointType> _constraintLipschitzConstants = std::vector<PointType>{})
        : IBaseConstrainedOptProblem<SearchAreaType, PointType>(_blockName, _problemNumber, _numberConstraints,
          _area, _optimalPoints, _optimalValue, _objectiveLipschitzConstant, _constraintLipschitzConstants),
          opt::IObjOptProblem<std::function<double(PointType, size_t)>>(_object) {};

    double computeObjectiveFunction(const PointType &x) const override { return this->object(x, this->numberConstraints); };
    double computeConstraintFunction(const PointType &x, size_t index) const override { return this->object(x, index); };
};

using OneDimensionalConstrainedOptProblem = ConstrainedOptProblem<opt::OneDimensionalSearchArea, double>;
using MultiDimensionalConstrainedOptProblem = ConstrainedOptProblem<opt::MultiDimensionalSearchArea, std::vector<double>>;

#endif // _CONSTRAINED_OPT_PROBLEM_H_
