#ifndef _ONE_DIMENSIONAL_FAMILY_OPT_PROBLEMS_GCGEN_H_
#define _ONE_DIMENSIONAL_FAMILY_OPT_PROBLEMS_GCGEN_H_

#include <opt_problems/base/IBaseFamilyOptProblems.h>
#include <general/classes/opt_problems/IObjOptProblem.h>
#include <general/structures/search_areas/OneDimensionalSearchArea.h>
#include <IOptProblemFamily.hpp>

class OneDimensionalFamilyOptProblemsGCGen :
    public IBaseFamilyOptProblems<opt::OneDimensionalSearchArea, double>, public opt::IObjOptProblem<IOptProblemFamily*>
{
public:
    OneDimensionalFamilyOptProblemsGCGen(IOptProblemFamily *_object = nullptr)
        : IBaseFamilyOptProblems<opt::OneDimensionalSearchArea, double>(),
          opt::IObjOptProblem<IOptProblemFamily*>(_object) {};

    void setObject(IOptProblemFamily* &_object) { object = _object; };
    void setFamilyName(std::string &_familyName) { familyName = _familyName; };

    size_t getFamilySize() const override { return object->GetFamilySize(); };

    opt::OneDimensionalSearchArea getSearchArea() const override {
        std::vector<double> A, B;
        (*object)[problemNumber]->GetBounds(A, B);
        return opt::OneDimensionalSearchArea(A[0], B[0]);
    };
    void getOptimalPoints(std::vector<double> &_optimalPoints) const override {
        _optimalPoints = (*object)[problemNumber]->GetOptimumPoint();
    };
    double getOptimalValue() const override {
        return (*object)[problemNumber]->GetOptimumValue();
    };
    double getObjectiveLipschitzConstant() const override { return (*object)[problemNumber]->GetLipschitzConstant(); };

    double computeObjectiveFunction(const double &x) const override {
        return (*object)[problemNumber]->ComputeFunction(std::vector<double>{ x });
    };
};

#endif // _ONE_DIMENSIONAL_FAMILY_OPT_PROBLEMS_GCGEN_H_
