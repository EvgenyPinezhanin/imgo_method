#ifndef ONE_DIMENSIONAL_PROBLEM_H_
#define ONE_DIMENSIONAL_PROBLEM_H_

#include <vector>
#include <functional>

#include <base_classes/opt_problems/IGeneralOptProblem.h>
#include <base_structures/search_areas/OneDimensionalSearchArea.h>

class OneDimensionalProblem:
    public opt::IGeneralOptProblem<std::function<double(double)>, opt::OneDimensionalSearchArea, double>
{
private:
    double lipschitzConstant;

public:
    OneDimensionalProblem(const std::function<double(double)> &_object = nullptr,
                          const opt::OneDimensionalSearchArea &_area = opt::OneDimensionalSearchArea(),
                          const std::vector<double> &_optimalPoints = std::vector<double>{},
                          double _optimalValue = 0.0,
                          double _lipschitzConstant = -1.0):
        opt::IGeneralOptProblem<std::function<double(double)>, opt::OneDimensionalSearchArea, double>(
            _object, _area, _optimalPoints, _optimalValue),
        lipschitzConstant(_lipschitzConstant)
    {};

    void setSearchArea(const opt::OneDimensionalSearchArea &_area) { area = _area; };
    void setOptimalPoints(const std::vector<double> &_optimalPoints) { optimalPoints = _optimalPoints; };
    void setOptimalValue(double _optimalValue) { optimalValue = _optimalValue; };

    void setLipschitzConstant(double _lipschitzConstant) { lipschitzConstant = _lipschitzConstant; };
    double getLipschitzConstant() const { return lipschitzConstant; };

    double computeObjFunction(const double &x) const override { return object(x); };
};

#endif // ONE_DIMENSIONAL_PROBLEM_H_
