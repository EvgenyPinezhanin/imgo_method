#ifndef ONE_DIMENSIONAL_PROBLEM_H_
#define ONE_DIMENSIONAL_PROBLEM_H_

#include <functional>

#include <base_classes/opt_problems/IGeneralOptProblem.h>
#include <base_classes/search_areas/OneDimensionalSearchArea.h>

using std::function;

class OneDimensionalProblem : public IGeneralOptProblem<function<double(double)>, OneDimensionalSearchArea, double> {
private:
    double lipschitzConstant;

public:
    OneDimensionalProblem(const function<double(double)> &_objFunction = nullptr,
                       const OneDimensionalSearchArea &_area = OneDimensionalSearchArea(0.0, 1.0),
                       vector<double> _optimalPoints = vector<double>{}, double _lipschitzConstant = -1.0)
                      : IGeneralOptProblem<function<double(double)>, OneDimensionalSearchArea, double>(_objFunction,
                      _area, _optimalPoints), lipschitzConstant(_lipschitzConstant) {};

    void setLipschitzConstant(double _lipschitzConstant) { lipschitzConstant = _lipschitzConstant; };
    double getLipschitzConstant() const { return lipschitzConstant; };

    double computeObjFunction(double x) const override { return objFunction(x); };
};

#endif // ONE_DIMENSIONAL_OPTIMIZATION_PROBLEM_H_
