#include <opt_problems/BaseFittingFamilyOptProblems.h>

#include <memory>

#include <MyMath.h>

double OneDimensionalSupportiveOptProblem::computeObjectiveFunction(const double &x) const {
    double result = 0.0;

    size_t omegaSize = omega.size();
    for (size_t i = 0; i < omegaSize; ++i) {
        result += coefficients[2 * i] * std::sin(omega[i] * x) +
                  coefficients[2 * i + 1] * std::cos(omega[i] * x);
    }

    return typeProblem == TypeProblem::MIN ? result : -result;
}

double BaseFittingFamilyOptProblems::uDerivative(const std::vector<double> &x) const {
    double result = 0.0;

    size_t numberX = x.size();
    for (size_t i = 0; i < numberX; ++i) {
        result += coefficients[2 * i] * x[i] * std::cos(x[i] * testPoints[numberTestPoints - 1].x[0]) -
                  coefficients[2 * i + 1] * x[i] * std::sin(x[i] * testPoints[numberTestPoints - 1].x[0]) ;
    }

    return result;
}

void BaseFittingFamilyOptProblems::calcCoefficients(const std::vector<double> &x) const {
    std::vector<std::vector<double>> A(numberCoefficients, std::vector<double>(numberCoefficients, 0));
    std::vector<double> B(numberCoefficients, 0);

    std::vector<std::function<double(double, double)>> functions{
        [] (double t, double x) -> double { return std::sin(x * t); },
        [] (double t, double x) -> double { return std::cos(x * t); }
    };

    mnk minimizer;

    for (size_t i = 0; i < numberCoefficients; ++i) {
        for (size_t j = 0; j < numberCoefficients; ++j) {
            for (size_t k = 0; k < numberTestPoints; ++k) {
                A[i][j] += functions[j % 2](testPoints[k].x[0], x[j / 2]) * functions[i % 2](testPoints[k].x[0], x[i / 2]);
            }
        }
        for (size_t j = 0; j < numberTestPoints; ++j) {
            B[i] += testPoints[j].x[1] * functions[i % 2](testPoints[j].x[0], x[i / 2]);
        }
    }

    minimizer.setA(A);
    minimizer.setB(B);
    minimizer.solve(coefficients);
}

bool BaseFittingFamilyOptProblems::equalPoints(
    const std::vector<double> &firstPoint, const std::vector<double> &secondPoint) const
{
    for (size_t i = 0; i < dimension; ++i) {
        if (firstPoint[i] != secondPoint[i]) {
            return false;
        }
    }
    return true;
}

void BaseFittingFamilyOptProblems::getOptimalPoints(std::vector<std::vector<double>> &_optimalPoints) const {
    if (!isSortX) {
        _optimalPoints = optimalPoints[problemNumber];
    } else {
        std::vector<size_t> indexes{ 0, 1, 2, 3 };
        _optimalPoints.resize(factorial(dimension));
        size_t index = 0;

        _optimalPoints[index] = optimalPoints[problemNumber][0];
        while (std::next_permutation(indexes.begin(), indexes.end())) {
            ++index;
            _optimalPoints[index] = std::vector<double>{ optimalPoints[problemNumber][0][indexes[0]],
                                                         optimalPoints[problemNumber][0][indexes[1]],
                                                         optimalPoints[problemNumber][0][indexes[2]],
                                                         optimalPoints[problemNumber][0][indexes[3]] };
        }
    }
}

double BaseFittingFamilyOptProblems::computeObjectiveFunction(const std::vector<double> &x) const {
    calcCoefficients(x);
    return -std::abs(uDerivative(x));
}

double BaseFittingFamilyOptProblems::computeConstraintFunction(const std::vector<double> &x, size_t index) const {
    double minValue, maxValue, value;

    if (index >= 0 && index < dimension - 1) {
        if (isSortX) {
            std::vector<double> sortX = x;
            std::sort(sortX.begin(), sortX.end());
            return sortX[index] * (1.0 + alpha) - sortX[index + 1] * (1.0 - alpha);
        }
        return x[index] * (1.0 + alpha) - x[index + 1] * (1.0 - alpha);
    } else if (index >= dimension - 1 && index < dimension + 2) {
        calcCoefficients(x);

        u.setOmega(x);
        u.setCoefficients(coefficients);
        u.setTypeProblem(OneDimensionalSupportiveOptProblem::TypeProblem::MIN);

        return std::abs(u.computeObjectiveFunction(testPoints[numberWindowPoints + index - dimension + 1].x[0]) -
                        testPoints[numberWindowPoints + index - dimension + 1].x[1]) - delta;
    } else if (index == dimension + 2) {
        calcCoefficients(x);

        u.setOmega(x);
        u.setCoefficients(coefficients);
        
        u.setTypeProblem(OneDimensionalSupportiveOptProblem::TypeProblem::MIN);

        minValue = u.computeObjectiveFunction(leftBoundWindow);
        for (double i = leftBoundWindow; i <= rightBoundWindow; i += 0.0005) {
            value = u.computeObjectiveFunction(i);
            if (value < minValue) {
                minValue = value;
            }
        }

        maxValue = u.computeObjectiveFunction(leftBoundWindow);
        for (double i = leftBoundWindow; i <= rightBoundWindow; i += 0.0005) {
            value = u.computeObjectiveFunction(i);
            if (value > maxValue) {
                maxValue = value;
            }
        }

        // using Result = typename GsaMethod<OneDimensionalSupportiveOptProblem>::Result;
        // std::unique_ptr<Result> result(static_cast<Result*>(gsa.createResult()));
        // 
        // u.setTypeProblem(OneDimensionalSupportiveOptProblem::TypeProblem::MIN);
        // gsa.setProblem(u);
        // gsa.solve(*result);
        // minValue = result->value;
        // 
        // u.setTypeProblem(OneDimensionalSupportiveOptProblem::TypeProblem::MAX);
        // gsa.setProblem(u);
        // gsa.solve(*result);
        // maxValue = -result->value;

        return maxValue - minValue - 2 * delta;
    } else if (index == dimension + 3) {
        calcCoefficients(x);
        return -std::abs(uDerivative(x));
    }
    return std::numeric_limits<double>::quiet_NaN();
}
