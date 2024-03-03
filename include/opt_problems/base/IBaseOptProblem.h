#ifndef _I_BASE_OPT_PROBLEM_H_
#define _I_BASE_OPT_PROBLEM_H_

#include <vector>

#include <general/classes/opt_problems/IGeneralOptProblem.h>

template <typename SearchAreaType, typename PointType>
class IBaseOptProblem : public opt::IGeneralOptProblem<SearchAreaType, PointType> {
protected:
    std::string blockName;
    size_t problemNumber;

    std::vector<PointType> optimalPoints;
    double optimalValue;

    double objectiveLipschitzConstant;

public:
    IBaseOptProblem(const std::string _blockName = "", size_t _problemNumber = 0,
                    const SearchAreaType &_area = SearchAreaType(),
                    const std::vector<PointType> &_optimalPoints = std::vector<PointType>{},
                    double _optimalValue = 0.0, double _objectiveLipschitzConstant = 0.0)
        : opt::IGeneralOptProblem<SearchAreaType, PointType>(_area), blockName(_blockName),
          problemNumber(_problemNumber), optimalPoints(_optimalPoints), optimalValue(_optimalValue),
          objectiveLipschitzConstant(_objectiveLipschitzConstant) {};

    void getBlockName(std::string &_blockName) const { _blockName = blockName; };
    size_t getProblemNumber() const { return problemNumber; };

    void getOptimalPoints(std::vector<PointType> &_optimalPoints) const override { _optimalPoints = optimalPoints; };
    double getOptimalValue() const override { return optimalValue; };

    virtual double getObjectiveLipschitzConstant() const { return objectiveLipschitzConstant; };
};

#endif // _I_BASE_OPT_PROBLEM_H_
