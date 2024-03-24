#ifndef _FITTING_TEST_FAMILY_OPT_PROBLEMS_H_
#define _FITTING_TEST_FAMILY_OPT_PROBLEMS_H_

#include <opt_problems/BaseFittingFamilyOptProblems.h>
#include <general/structures/search_areas/MultiDimensionalSearchArea.h>

class FittingFamilyOptProblems : public BaseFittingFamilyOptProblems {
public:
    FittingFamilyOptProblems(
        const GsaMethod<OneDimensionalSupportiveOptProblem> &_gsa = GsaMethod<OneDimensionalSupportiveOptProblem>());
};

#endif // _FITTING_TEST_FAMILY_OPT_PROBLEMS_H_
