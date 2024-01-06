#ifndef ONE_DIMENSIONAL_CONSTRAINED_PROBLEM_H_
#define ONE_DIMENSIONAL_CONSTRAINED_PROBLEM_H_

#include <vector>
#include <functional>

#include <base_classes/opt_problems/IGeneralConstrainedOptProblem.h>
#include <base_structures/search_areas/OneDimensionalSearchArea.h>

class OneDimensionalConstrainedProblem:
    public opt::IGeneralConstrainedOptProblem<std::function<double(double, int)>,
                                              opt::OneDimensionalSearchArea, double>
{
private:
    std::string blockName;
    size_t problemNumber;

    std::vector<double> lipschitzConstants;

public:
    OneDimensionalConstrainedProblem(const std::function<double(double, int)> &_object = nullptr, const std::string &_blockName = "",
                                     size_t _problemNumber = 0,
                                     const opt::OneDimensionalSearchArea &_area = opt::OneDimensionalSearchArea(),
                                     int _numberConstraints = 0,
                                     const std::vector<double> &_optimalPoints = std::vector<double>{},
                                     double _optimalValue = 0.0,
                                     const std::vector<double> _lipschitzConstants = std::vector<double>{})
        : opt::IGeneralConstrainedOptProblem<std::function<double(double, int)>, opt::OneDimensionalSearchArea, double>(
          _object, _area, _numberConstraints, _optimalPoints, _optimalValue), blockName(_blockName), problemNumber(_problemNumber),
          lipschitzConstants(_lipschitzConstants) {};

    void setSearchArea(const opt::OneDimensionalSearchArea &_area) { area = _area; };
    void setOptimalPoints(const std::vector<double> &_optimalPoints) { optimalPoints = _optimalPoints; };
    void setOptimalValue(double _optimalValue) { optimalValue = _optimalValue; };
    void setNumberConstraints(size_t _numberConstraints) { numberConstraints = _numberConstraints; };

    void setBlockName(const std::string &_blockName) { blockName = _blockName; };
    void getBlockName(std::string &_blockName) const { _blockName = blockName; };

    void setProblemNumber(size_t _problemNumber) { problemNumber = _problemNumber; };
    size_t getProblemNumber() const { return problemNumber; };

    void setLipschitzConstants(const std::vector<double> &_lipschitzConstants) { lipschitzConstants = _lipschitzConstants; };
    void getLipschitzConstants(std::vector<double> &_lipschitzConstants) const { _lipschitzConstants = lipschitzConstants; };

    double computeObjFunction(const double &x) const override { return object(x, numberConstraints); };
    double computeConstraint(const double &x, int index) const override { return object(x, index); };
};

#endif // ONE_DIMENSIONAL_CONSTRAINED_PROBLEM_H_
