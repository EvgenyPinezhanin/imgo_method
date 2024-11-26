#ifndef _BASE_FITTING_FAMILY_OPT_PROBLEMS_H_
#define _BASE_FITTING_FAMILY_OPT_PROBLEMS_H_

#include <vector>
#include <functional>
#include <memory>

#include <general/classes/opt_methods/IGeneralNumericalOptMethod.h>
#include <opt_problems/base/IBaseOptProblem.h>
#include <general/structures/search_areas/OneDimensionalSearchArea.h>
#include <opt_problems/base/IBaseConstrainedFamilyOptProblems.h>
#include <general/structures/search_areas/MultiDimensionalSearchArea.h>
#include <MyMath.h>

class OneDimensionalSupportiveOptProblem : public IBaseOptProblem<opt::OneDimensionalSearchArea, double> {
public:
    enum class TypeProblem { MIN, MAX };

protected:
    std::vector<double> coefficients, omega;
    TypeProblem typeProblem;

public:
    OneDimensionalSupportiveOptProblem(const opt::OneDimensionalSearchArea &_area = opt::OneDimensionalSearchArea(),
                                       TypeProblem _typeProblem = TypeProblem::MIN)
        : IBaseOptProblem<opt::OneDimensionalSearchArea, double>("", 0, _area, std::vector<double>{}, 0.0, 0.0),
          coefficients(), omega(), typeProblem(_typeProblem) {};

    void setCoefficients(const std::vector<double> &_coefficients) { coefficients = _coefficients; };
    void getCoefficients(std::vector<double> &_coefficients) const { _coefficients = coefficients; };

    void setOmega(const std::vector<double> &_omega) { omega = _omega; };
    void getOmega(std::vector<double> &_omega) const { _omega = omega; };

    void setTypeProblem(TypeProblem _typeProblem) { typeProblem = _typeProblem; };

    double computeObjectiveFunction(const double &x) const override;
};

template<typename OptMethod>
class BaseFittingFamilyOptProblems :
    public IBaseConstrainedFamilyOptProblems<opt::MultiDimensionalSearchArea, std::vector<double>>
{
protected:
    size_t dimension;

    double alpha, delta;
    double leftBoundWideWindow, rightBoundWideWindow;
    std::vector<double> firstShortWindowValues, secondShortWindowValues;
    size_t numberWideWindowPoints, numberTestPoints;
    mutable std::vector<point> testPoints;

    size_t numberCoefficients;
    mutable std::vector<double> coefficients;

    mutable OptMethod optMethod;
    mutable OneDimensionalSupportiveOptProblem u;

    bool isSortX;

    double uDerivative(const std::vector<double> &x) const;

    void calcCoefficients(const std::vector<double> &x) const;

    bool equalPoints(const std::vector<double> &firstPoint, const std::vector<double> &secondPoint) const;

public:
    BaseFittingFamilyOptProblems(
        size_t _familySize = 0, size_t _dimension = 1,
        const opt::MultiDimensionalSearchArea &_area = opt::MultiDimensionalSearchArea(),
        double _alpha = 0.0, double _delta = 0.0,
        double _leftBoundWideWindow = 0.0, double _rightBoundWideWindow = 0.0,
        const std::vector<double> &insideWideWindowPoints = std::vector<double>{},
        double firstShortWindowPoint = 0.0, const std::vector<double> &_firstShortWindowValues = std::vector<double>{},
        double secondShortWindowPoint = 0.0, const std::vector<double> &_secondShortWindowValues = std::vector<double>{},
        double lastShortWindowPoint = 0.0,
        const OptMethod  &_optMethod = OptMethod(),
        bool _isSortX = false,
        const std::vector<std::vector<std::vector<double>>> &_optimalPoints = std::vector<std::vector<std::vector<double>>>{},
        const std::vector<double> &_optimalValue = std::vector<double>{},
        const std::vector<double> &_objectiveLipschitzConstant = std::vector<double>{},
        const std::vector<std::vector<double>> &_constraintLipschitzConstants = std::vector<std::vector<double>>{})
        : IBaseConstrainedFamilyOptProblems<opt::MultiDimensionalSearchArea, std::vector<double>>("Fitting Family",
          _familySize, _dimension + 3, _area, _optimalPoints, _optimalValue, _objectiveLipschitzConstant,
          _constraintLipschitzConstants),
          dimension(_dimension),
          alpha(_alpha), delta(_delta),
          leftBoundWideWindow(_leftBoundWideWindow), rightBoundWideWindow(_rightBoundWideWindow),
          firstShortWindowValues(_firstShortWindowValues), secondShortWindowValues(_secondShortWindowValues),
          numberWideWindowPoints(insideWideWindowPoints.size() + 2),
          numberTestPoints(numberWideWindowPoints + 3),
          testPoints(numberTestPoints),
          numberCoefficients(2 * dimension),
          coefficients(numberCoefficients),
          optMethod(_optMethod),
          u(opt::OneDimensionalSearchArea(leftBoundWideWindow, rightBoundWideWindow)),
          isSortX(_isSortX)
    {
        testPoints[0] = point(leftBoundWideWindow, 0.0);
        for (size_t i = 0; i < numberWideWindowPoints - 2; ++i) {
            testPoints[i + 1] = point(insideWideWindowPoints[i], 0.0);
        }
        testPoints[numberTestPoints - 4] = point(rightBoundWideWindow, 0.0);

        testPoints[numberTestPoints - 3] = point(firstShortWindowPoint, firstShortWindowValues[0]);
        testPoints[numberTestPoints - 2] = point(secondShortWindowPoint, secondShortWindowValues[0]);

        testPoints[numberTestPoints - 1] = point(lastShortWindowPoint, 0.0);
    };

    void setProblemNumber(size_t _problemNumber) const override {
        problemNumber = _problemNumber;
        testPoints[numberTestPoints - 3].x[1] = firstShortWindowValues[problemNumber];
        testPoints[numberTestPoints - 2].x[1] = secondShortWindowValues[problemNumber];
    };

    size_t getDimension() const { return dimension; };

    double getAlpha() const { return alpha; };
    double getDelta() const { return delta; };

    double getLeftBoundWideWindow() const { return leftBoundWideWindow; };
    double getRightBoundWideWindow() const { return rightBoundWideWindow; };
    size_t getNumberWideWindowPoints() const { return numberWideWindowPoints; };

    void getFirstShortWindowValues(std::vector<double> &_firstValues) const { _firstValues = firstShortWindowValues; };
    void getSecondShortWindowValues(std::vector<double> &_secondValues) const { _secondValues = secondShortWindowValues; };

    size_t getNumberTestPoints() const { return numberTestPoints; };
    void getTestPoints(std::vector<point> &_testPoints) const { _testPoints = testPoints; };

    size_t getNumberCoefficients() const { return numberCoefficients; };
    void getCoefficients(std::vector<double> &_coefficients) const { _coefficients = coefficients; };

    void getOptMethod(OptMethod &_optMethod) const { _optMethod = optMethod; };
    void setOptMethod(const OptMethod &_optMethod) {
        optMethod = _optMethod;
    };

    void setIsSortX(bool _isSortX) { isSortX = _isSortX; };

    void getOptimalPoints(std::vector<std::vector<double>> &_optimalPoints) const override;

    double computeObjectiveFunction(const std::vector<double> &x) const override;
    double computeConstraintFunction(const std::vector<double> &x, size_t index) const override;
};

template<typename OptMethod>
double BaseFittingFamilyOptProblems<OptMethod>::uDerivative(const std::vector<double> &x) const {
    double result = 0.0;

    size_t numberX = x.size();
    for (size_t i = 0; i < numberX; ++i) {
        result += coefficients[2 * i] * x[i] * std::cos(x[i] * testPoints[numberTestPoints - 1].x[0]) -
                  coefficients[2 * i + 1] * x[i] * std::sin(x[i] * testPoints[numberTestPoints - 1].x[0]) ;
    }

    return result;
}

template<typename OptMethod>
void BaseFittingFamilyOptProblems<OptMethod>::calcCoefficients(const std::vector<double> &x) const {
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

    for (size_t i = 0; i < numberCoefficients; ++i) {
        for (size_t j = 0; j < numberCoefficients; ++j) {
            std::cout << std::setprecision(10) << A[i][j] << " ";
        }
        std::cout << std::setprecision(10)  << B[i] << "\n";
    }

    minimizer.setA(A);
    minimizer.setB(B);
    minimizer.solve(coefficients);
}

template<typename OptMethod>
bool BaseFittingFamilyOptProblems<OptMethod>::equalPoints(
    const std::vector<double> &firstPoint, const std::vector<double> &secondPoint) const
{
    for (size_t i = 0; i < dimension; ++i) {
        if (firstPoint[i] != secondPoint[i]) {
            return false;
        }
    }
    return true;
}

template<typename OptMethod>
void BaseFittingFamilyOptProblems<OptMethod>::getOptimalPoints(std::vector<std::vector<double>> &_optimalPoints) const {
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

template<typename OptMethod>
double BaseFittingFamilyOptProblems<OptMethod>::computeObjectiveFunction(const std::vector<double> &x) const {
    calcCoefficients(x);
    return -std::abs(uDerivative(x));
}

template<typename OptMethod>
double BaseFittingFamilyOptProblems<OptMethod>::computeConstraintFunction(const std::vector<double> &x, size_t index) const {
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

        return std::abs(u.computeObjectiveFunction(testPoints[numberWideWindowPoints + index - dimension + 1].x[0]) -
                        testPoints[numberWideWindowPoints + index - dimension + 1].x[1]) - delta;
    } else if (index == dimension + 2) {
        calcCoefficients(x);

        u.setOmega(x);
        u.setCoefficients(coefficients);

        using Result = typename OptMethod::Result;
        std::unique_ptr<Result> result(static_cast<Result*>(optMethod.createResult()));
        
        u.setTypeProblem(OneDimensionalSupportiveOptProblem::TypeProblem::MIN);
        optMethod.setProblem(u);
        optMethod.solve(*result);
        minValue = result->value;
        
        u.setTypeProblem(OneDimensionalSupportiveOptProblem::TypeProblem::MAX);
        optMethod.setProblem(u);
        optMethod.solve(*result);
        maxValue = -result->value;

        // std::cout << maxValue - minValue - 2 * delta << "\n";

        return maxValue - minValue - 2 * delta;
    } else if (index == dimension + 3) {
        calcCoefficients(x);
        return -std::abs(uDerivative(x));
    }
    return std::numeric_limits<double>::quiet_NaN();
}

#endif // _BASE_FITTING_FAMILY_OPT_PROBLEMS_H_
