#ifndef _OPT_PROBLEM_H_
#define _OPT_PROBLEM_H_

#include <string>
#include <vector>
#include <functional>

#include <opt_problems/base/IBaseOptProblem.h>
#include <general/classes/opt_problems/IObjOptProblem.h>
#include <general/structures/search_areas/OneDimensionalSearchArea.h>
#include <general/structures/search_areas/MultiDimensionalSearchArea.h>

template <typename SearchAreaType, typename PointType>
class OptProblem :
    public IBaseOptProblem<SearchAreaType, PointType>, public opt::IObjOptProblem<std::function<double(PointType)>>
{
public:
    OptProblem(const std::function<double(PointType)> &_object = nullptr,
               const std::string &_blockName = "", size_t _problemNumber = 0,
               const SearchAreaType &_area = SearchAreaType(),
               const std::vector<PointType> &_optimalPoints = std::vector<PointType>{},
               double _optimalValue = 0.0, double _objectiveLipschitzConstant = 0.0)
        : IBaseOptProblem<SearchAreaType, PointType>(_blockName, _problemNumber, _area,
          _optimalPoints, _optimalValue, _objectiveLipschitzConstant),
          opt::IObjOptProblem<std::function<double(PointType)>>(_object) {};

    double computeObjectiveFunction(const PointType &x) const override { return this->object(x); };
};

using OneDimensionalOptProblem = OptProblem<opt::OneDimensionalSearchArea, double>;
using MultiDimensionalOptProblem = OptProblem<opt::MultiDimensionalSearchArea, std::vector<double>>;

#endif // _OPT_PROBLEM_H_
