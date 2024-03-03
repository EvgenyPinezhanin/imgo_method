#ifndef _I_BASE_CONSTRAINED_OPT_PROBLEM_H_
#define _I_BASE_CONSTRAINED_OPT_PROBLEM_H_

#include <general/classes/opt_problems/IConstrainedOptProblem.h>
#include <opt_problems/base/IBaseOptProblem.h>

template <typename SearchAreaType, typename PointType>
class IBaseConstrainedOptProblem :
    public IBaseOptProblem<SearchAreaType, PointType>, public opt::IConstrainedOptProblem<PointType>
{
protected:
    std::vector<double> constraintLipschitzConstants;

public:
    IBaseConstrainedOptProblem(const std::string _blockName = "", size_t _problemNumber = 0,
                               size_t _numberConstraints = 0, const SearchAreaType &_area = SearchAreaType(),
                               const std::vector<PointType> &_optimalPoints = std::vector<PointType>{},
                               double _optimalValue = 0.0, double _objectiveLipschitzConstant = 0.0,
                               const std::vector<double> &_constraintLipschitzConstants = std::vector<double>{})
        : IBaseOptProblem<SearchAreaType, PointType>(_blockName, _problemNumber, _area, _optimalPoints,
          _optimalValue, _objectiveLipschitzConstant),
          opt::IConstrainedOptProblem<PointType>(_numberConstraints),
          constraintLipschitzConstants(_constraintLipschitzConstants) {};

    virtual void getConstraintLipschitzConstants(std::vector<double> _constraintLipschitzConstants) const {
        _constraintLipschitzConstants = constraintLipschitzConstants;
    };
};

#endif // _I_BASE_CONSTRAINED_OPT_PROBLEM_H_
