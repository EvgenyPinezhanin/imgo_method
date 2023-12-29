#ifndef MULTI_DIMENSIONAL_CONSTRAINED_PROBLEM_H_
#define MULTI_DIMENSIONAL_CONSTRAINED_PROBLEM_H_

#include <vector>
#include <functional>

#include <base_classes/opt_problems/IGeneralConstrainedOptProblem.h>
#include <base_structures/search_areas/MultiDimensionalSearchArea.h>

class MultiDimensionalConstrainedProblem:
    public opt::IGeneralConstrainedOptProblem<std::function<double(std::vector<double>, int)>,
                                              opt::MultiDimensionalSearchArea, std::vector<double>>
{
private:
    std::vector<double> lipschitzConstants;

public:
    MultiDimensionalConstrainedProblem(const std::function<double(std::vector<double>, int)> &_object = nullptr,
                                       const opt::MultiDimensionalSearchArea &_area = opt::MultiDimensionalSearchArea(),
                                       int _numberConstraints = 0,
                                       const std::vector<std::vector<double>> &_optimalPoints = std::vector<std::vector<double>>{},
                                       double _optimalValue = 0.0,
                                       const std::vector<double> _lipschitzConstants = std::vector<double>{})
        : opt::IGeneralConstrainedOptProblem<std::function<double(std::vector<double>, int)>,
                                             opt::MultiDimensionalSearchArea, std::vector<double>>(
          _object, _area, _numberConstraints, _optimalPoints, _optimalValue), lipschitzConstants(_lipschitzConstants) {};

    void setSearchArea(const opt::MultiDimensionalSearchArea &_area) { area = _area; };
    void setOptimalPoints(const std::vector<std::vector<double>> &_optimalPoints) { optimalPoints = _optimalPoints; };
    void setOptimalValue(double _optimalValue) { optimalValue = _optimalValue; };

    void setLipschitzConstants(const std::vector<double> &_lipschitzConstants) { lipschitzConstants = _lipschitzConstants; };
    void getLipschitzConstants(std::vector<double> &_lipschitzConstants) const { _lipschitzConstants = lipschitzConstants; };

    double computeObjFunction(const std::vector<double> &x) const override { return object(x, numberConstraints); };
    double computeConstraint(const std::vector<double> &x, int index) const override { return object(x, index); };
};

#endif // MULTI_DIMENSIONAL_CONSTRAINED_PROBLEM_H_
