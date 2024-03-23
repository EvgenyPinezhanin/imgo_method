#ifndef _I_BASE_CONSTRAINED_FAMILY_OPT_PROBLEMS_H_
#define _I_BASE_CONSTRAINED_FAMILY_OPT_PROBLEMS_H_

#include <vector>
#include <string>

#include <opt_problems/base/IBaseFamilyOptProblems.h>
#include <general/classes/opt_problems/IConstrainedOptProblem.h>

template <typename SearchAreaType, typename PointType>
class IBaseConstrainedFamilyOptProblems :
    public IBaseFamilyOptProblems<SearchAreaType, PointType>, public opt::IConstrainedOptProblem<PointType>
{
protected:
    std::vector<std::vector<double>> constraintLipschitzConstants;

public:
    IBaseConstrainedFamilyOptProblems(
        const std::string _familyName = "", size_t _familySize = 0, size_t _numberConstraints = 0,
        const SearchAreaType &_area = SearchAreaType(),
        const std::vector<std::vector<PointType>> &_optimalPoints = std::vector<std::vector<PointType>>{},
        const std::vector<double> &_optimalValue = std::vector<double>{},
        const std::vector<double> &_objectiveLipschitzConstant = std::vector<double>{},
        const std::vector<std::vector<double>> &_constraintLipschitzConstants = std::vector<std::vector<double>>{})
        : IBaseFamilyOptProblems<SearchAreaType, PointType>(_familyName, _familySize, _area, _optimalPoints,
          _optimalValue, _objectiveLipschitzConstant),
          opt::IConstrainedOptProblem<PointType>(_numberConstraints),
          constraintLipschitzConstants(_constraintLipschitzConstants) {};

    virtual void getConstraintLipschitzConstants(std::vector<double> _constraintLipschitzConstants) const {
        _constraintLipschitzConstants = constraintLipschitzConstants[this->problemNumber];
    };
};

#endif // _I_BASE_CONSTRAINED_FAMILY_OPT_PROBLEMS_H_
