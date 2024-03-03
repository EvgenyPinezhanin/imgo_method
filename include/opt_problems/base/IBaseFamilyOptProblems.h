#ifndef _I_BASE_FAMILY_OPT_PROBLEMS_H_
#define _I_BASE_FAMILY_OPT_PROBLEMS_H_

#include <vector>

#include <general/classes/opt_problems/IGeneralOptProblem.h>

template <typename SearchAreaType, typename PointType>
class IBaseFamilyOptProblems : public opt::IGeneralOptProblem<SearchAreaType, PointType> {
protected:
    std::string familyName;
    size_t familySize;
    mutable size_t problemNumber;

    std::vector<std::vector<PointType>> optimalPoints;
    std::vector<double> optimalValue;

    std::vector<double> objectiveLipschitzConstant;

public:
    IBaseFamilyOptProblems(const std::string _familyName = "", size_t _familySize = 0,
                           const SearchAreaType &_area = SearchAreaType(),
                           const std::vector<std::vector<PointType>> &_optimalPoints = std::vector<std::vector<PointType>>{},
                           const std::vector<double> &_optimalValue = std::vector<double>{},
                           const std::vector<double> &_objectiveLipschitzConstant = std::vector<double>{})
        : opt::IGeneralOptProblem<SearchAreaType, PointType>(_area), familyName(_familyName),
          familySize(_familySize), problemNumber(0), optimalPoints(_optimalPoints),
          optimalValue(_optimalValue), objectiveLipschitzConstant(_objectiveLipschitzConstant) {};

    void getFamilyName(std::string &_familyName) const { _familyName = familyName; };
    virtual size_t getFamilySize() const { return familySize; };

    void setProblemNumber(size_t _problemNumber) const { problemNumber = _problemNumber; };
    size_t getProblemNumber() const { return problemNumber; };

    virtual void getOptimalPoints(std::vector<PointType> &_optimalPoints) const override {
        _optimalPoints = optimalPoints[problemNumber];
    };
    double getOptimalValue() const override { return optimalValue[problemNumber]; };

    virtual double getObjectiveLipschitzConstant() const { return objectiveLipschitzConstant[problemNumber]; };
};

#endif // _I_BASE_FAMILY_OPT_PROBLEMS_H_
