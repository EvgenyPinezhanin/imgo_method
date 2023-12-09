#ifndef ONE_DIM_CONSTRAINED_PROBLEM_H_
#define ONE_DIM_CONSTRAINED_PROBLEM_H_

#include <vector>
#include <functional>

#include <base_classes/opt_problems/IGeneralConstrainedOptProblem.h>
#include <base_structures/search_areas/OneDimensionalSearchArea.h>

class OneDimConstrainedProblem:
    public opt::IGeneralConstrainedOptProblem<std::function<double(double, int)>,
                                              opt::OneDimensionalSearchArea, double>
{
private:
    std::vector<double> lipschitzConstants;

public:
    OneDimConstrainedProblem(const std::function<double(double, int)> &_object = nullptr,
                             const opt::OneDimensionalSearchArea &_area = opt::OneDimensionalSearchArea(),
                             int _numberConstraints = 0,
                             const std::vector<double> &_optimalPoints = std::vector<double>{},
                             double _optimalValue = 0.0,
                             const std::vector<double> _lipschitzConstants = std::vector<double>{})
        : opt::IGeneralConstrainedOptProblem<std::function<double(double, int)>, opt::OneDimensionalSearchArea, double>(
          _object, _area, _numberConstraints, _optimalPoints, _optimalValue), lipschitzConstants(_lipschitzConstants) {};

    void setSearchArea(const opt::OneDimensionalSearchArea &_area) { area = _area; };
    void setOptimalPoints(const std::vector<double> &_optimalPoints) { optimalPoints = _optimalPoints; };
    void setOptimalValue(double _optimalValue) { optimalValue = _optimalValue; };

    void setLipschitzConstants(const std::vector<double> &_lipschitzConstants) { lipschitzConstants = _lipschitzConstants; };
    void getLipschitzConstants(std::vector<double> &_lipschitzConstants) const { _lipschitzConstants = lipschitzConstants; };

    double computeObjFunction(const double &x) const override { return object(x, numberConstraints); };
    double computeConstraint(const double &x, int index) const override { return object(x, index); };
};

#endif // ONE_DIM_CONSTRAINED_PROBLEM_H_
