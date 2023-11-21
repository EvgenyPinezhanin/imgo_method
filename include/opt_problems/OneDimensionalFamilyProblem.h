#ifndef ONE_DIMENSIONAL_FAMILY_PROBLEM_H_
#define ONE_DIMENSIONAL_FAMILY_PROBLEM_H_

#include <base_classes/opt_problems/IGeneralOptProblem.h>
#include <base_structures/search_areas/OneDimensionalSearchArea.h>
#include <IOptProblemFamily.hpp>

class OneDimensionalFamilyProblem:
    public opt::IGeneralOptProblem<IOptProblemFamily>, opt::OneDimensionalSearchArea, double>
{
private:
    double lipschitzConstant;
    int currentProblem;

public:
    OneDimensionalFamilyProblem(const function<double(double)> &_objFunction = nullptr,
                          const opt::OneDimensionalSearchArea &_area = opt::OneDimensionalSearchArea(0.0, 1.0),
                          vector<double> _optimalPoints = vector<double>{}, double _optimalValue = 0.0,
                          double _lipschitzConstant = -1.0)
                         : opt::IGeneralOptProblem<function<double(double)>, opt::OneDimensionalSearchArea, double>(_objFunction,
                         _area, _optimalPoints, _optimalValue), lipschitzConstant(_lipschitzConstant) {};

    void setLipschitzConstant(double _lipschitzConstant) { lipschitzConstant = _lipschitzConstant; };
    double getLipschitzConstant() const { return lipschitzConstant; };

    double computeObjFunction(const double &x) const override { return objFunction(x); };
};

#endif // ONE_DIMENSIONAL_FAMILY_PROBLEM_H_
