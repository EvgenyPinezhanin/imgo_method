#ifndef MULTI_DIMENSIONAL_CONSTRAINED_FITTING_FAMILY_PROBLEM_H_
#define MULTI_DIMENSIONAL_CONSTRAINED_FITTING_FAMILY_PROBLEM_H_

#include <vector>

#include <base_classes/opt_problems/IGeneralConstrainedOptProblem.h>
#include <base_structures/search_areas/MultiDimensionalSearchArea.h>

/* class OneDimensionalFamilyProblemsGCGen:
    public opt::IGeneralConstrainedOptProblem<IOptProblemFamily*, opt::MultiDimensionalSearchArea, std::vector<double>>
{
private:
    mutable size_t currentProblem;

public:
    OneDimensionalFamilyProblemsGCGen(IOptProblemFamily *_object = nullptr)
        : opt::IGeneralOptProblem<IOptProblemFamily*, opt::OneDimensionalSearchArea, double>(_object),
          currentProblem(0) {};

    void setCurrentProblem(size_t _currentProblem) const { currentProblem = _currentProblem; };
    size_t getCurrentProblem() const { return currentProblem; }

    size_t getFamilySize() const { return object->GetFamilySize(); };

    opt::OneDimensionalSearchArea getSearchArea() const override {
        std::vector<double> A, B;
        (*object)[currentProblem]->GetBounds(A, B);
        return opt::OneDimensionalSearchArea(A[0], B[0]);
    };
    void getOptimalPoints(std::vector<double> &_optimalPoints) const override {
        _optimalPoints = (*object)[currentProblem]->GetOptimumPoint();
    };
    double getOptimalValue() const override {
        return (*object)[currentProblem]->GetOptimumValue();
    };
    double getLipschitzConstant() const { return (*object)[currentProblem]->GetLipschitzConstant(); };

    double computeObjFunction(const double &x) const override {
        return (*object)[currentProblem]->ComputeFunction(std::vector<double>{x});
    };
}; */

#endif // MULTI_DIMENSIONAL_CONSTRAINED_FITTING_FAMILY_PROBLEM_H_
