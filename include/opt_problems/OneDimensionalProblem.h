#ifndef ONE_DIMENSIONAL_PROBLEM_H_
#define ONE_DIMENSIONAL_PROBLEM_H_

#include <string>
#include <vector>
#include <functional>

#include <base_classes/opt_problems/IGeneralOptProblem.h>
#include <base_structures/search_areas/OneDimensionalSearchArea.h>

class OneDimensionalProblem:
    public opt::IGeneralOptProblem<std::function<double(double)>, opt::OneDimensionalSearchArea, double>
{
private:
    std::string blockName;
    size_t problemNumber;

    double lipschitzConstant;

public:
    OneDimensionalProblem(const std::function<double(double)> &_object = nullptr, const std::string &_blockName = "",
                          size_t _problemNumber = 0,
                          const opt::OneDimensionalSearchArea &_area = opt::OneDimensionalSearchArea(),
                          const std::vector<double> &_optimalPoints = std::vector<double>{},
                          double _optimalValue = 0.0, double _lipschitzConstant = -1.0)
        : opt::IGeneralOptProblem<std::function<double(double)>, opt::OneDimensionalSearchArea, double>(
          _object, _area, _optimalPoints, _optimalValue), blockName(_blockName), problemNumber(_problemNumber),
          lipschitzConstant(_lipschitzConstant) {};

    void setSearchArea(const opt::OneDimensionalSearchArea &_area) { area = _area; };
    void setOptimalPoints(const std::vector<double> &_optimalPoints) { optimalPoints = _optimalPoints; };
    void setOptimalValue(double _optimalValue) { optimalValue = _optimalValue; };

    void setBlockName(const std::string &_blockName) { blockName = _blockName; };
    void getBlockName(std::string &_blockName) const { _blockName = blockName; };

    void setProblemNumber(size_t _problemNumber) { problemNumber = _problemNumber; };
    size_t getProblemNumber() const { return problemNumber; };

    void setLipschitzConstant(double _lipschitzConstant) { lipschitzConstant = _lipschitzConstant; };
    double getLipschitzConstant() const { return lipschitzConstant; };

    double computeObjFunction(const double &x) const override { return object(x); };
};

#endif // ONE_DIMENSIONAL_PROBLEM_H_
