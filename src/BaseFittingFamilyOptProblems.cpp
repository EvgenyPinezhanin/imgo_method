#include <opt_problems/BaseFittingFamilyOptProblems.h>

double OneDimensionalSupportiveOptProblem::computeObjectiveFunction(const double &x) const {
    double result = 0.0;

    size_t omegaSize = omega.size();
    for (size_t i = 0; i < omegaSize; ++i) {
        result += coefficients[2 * i] * std::sin(omega[i] * x) +
                  coefficients[2 * i + 1] * std::cos(omega[i] * x);
    }

    return typeProblem == TypeProblem::MIN ? result : -result;
}
