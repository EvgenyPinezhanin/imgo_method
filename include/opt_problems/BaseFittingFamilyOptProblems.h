#ifndef _BASE_FITTING_FAMILY_OPT_PROBLEMS_H_
#define _BASE_FITTING_FAMILY_OPT_PROBLEMS_H_

#include <vector>
#include <functional>

#include <opt_problems/base/IBaseOptProblem.h>
#include <general/structures/search_areas/OneDimensionalSearchArea.h>
#include <opt_problems/base/IBaseConstrainedFamilyOptProblems.h>
#include <general/structures/search_areas/MultiDimensionalSearchArea.h>
#include <opt_methods/GsaMethod.h>

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

class BaseFittingFamilyOptProblems :
    public IBaseConstrainedFamilyOptProblems<opt::MultiDimensionalSearchArea, std::vector<double>>
{
protected:
    size_t dimension;

    double alpha, delta, leftBoundWindow, rightBoundWindow;
    std::vector<double> firstValues, secondValues;
    size_t numberWindowPoints, numberTestPoints;
    mutable std::vector<point> testPoints;

    size_t numberCoefficients;
    mutable std::vector<double> coefficients;

    // TODO: think about mutable
    mutable GsaMethod<OneDimensionalSupportiveOptProblem> gsa;
    mutable OneDimensionalSupportiveOptProblem u;

    bool isSortX;

    double uDerivative(const std::vector<double> &x) const;

    void calcCoefficients(const std::vector<double> &x) const;

    bool equalPoints(const std::vector<double> &firstPoint, const std::vector<double> &secondPoint) const;

public:
    BaseFittingFamilyOptProblems(
        size_t _familySize = 0, size_t _dimension = 1,
        const opt::MultiDimensionalSearchArea &_area = opt::MultiDimensionalSearchArea(),
        double _alpha = 0.0, double _delta = 0.0, double _leftBoundWindow = 0.0, double _rightBoundWindow = 0.0,
        double firstPoint = 0.0, const std::vector<double> &_firstValues = std::vector<double>{},
        double secondPoint = 0.0, const std::vector<double> &_secondValues = std::vector<double>{}, double lastPoint = 0.0,
        const std::vector<double> &windowPoints = std::vector<double>{},
        const GsaMethod<OneDimensionalSupportiveOptProblem> &_gsa = GsaMethod<OneDimensionalSupportiveOptProblem>(), bool _isSortX = false,
        const std::vector<std::vector<std::vector<double>>> &_optimalPoints = std::vector<std::vector<std::vector<double>>>{},
        const std::vector<double> &_optimalValue = std::vector<double>{},
        const std::vector<double> &_objectiveLipschitzConstant = std::vector<double>{},
        const std::vector<std::vector<double>> &_constraintLipschitzConstants = std::vector<std::vector<double>>{})
        : IBaseConstrainedFamilyOptProblems<opt::MultiDimensionalSearchArea, std::vector<double>>("Fitting Family",
          _familySize, _dimension + 3, _area, _optimalPoints, _optimalValue, _objectiveLipschitzConstant,
          _constraintLipschitzConstants), dimension(_dimension), alpha(_alpha), delta(_delta),
          leftBoundWindow(_leftBoundWindow), rightBoundWindow(_rightBoundWindow), firstValues(_firstValues),
          secondValues(_secondValues), numberWindowPoints(windowPoints.size() + 2), numberTestPoints(numberWindowPoints + 3),
          testPoints(numberTestPoints), numberCoefficients(2 * dimension), coefficients(numberCoefficients), gsa(_gsa),
          u(opt::OneDimensionalSearchArea(leftBoundWindow, rightBoundWindow)), isSortX(_isSortX)
    {
        testPoints[0] = point(leftBoundWindow, 0.0);
        for (size_t i = 0; i < numberWindowPoints - 2; ++i) {
            testPoints[i + 1] = point(windowPoints[i], 0.0);
        }
        testPoints[numberTestPoints - 4] = point(rightBoundWindow, 0.0);
        testPoints[numberTestPoints - 3] = point(firstPoint, firstValues[0]);
        testPoints[numberTestPoints - 2] = point(secondPoint, secondValues[0]);
        testPoints[numberTestPoints - 1] = point(lastPoint, 0.0);
    };

    void setProblemNumber(size_t _problemNumber) const override {
        problemNumber = _problemNumber;
        testPoints[numberTestPoints - 3].x[1] = firstValues[problemNumber];
        testPoints[numberTestPoints - 2].x[1] = secondValues[problemNumber];
    };

    size_t getDimension() const { return dimension; };

    double getAlpha() const { return alpha; };
    double getDelta() const { return delta; };

    double getLeftBoundWindow() const { return leftBoundWindow; };
    double getRightBoundWindow() const { return rightBoundWindow; };

    void getFirstValues(std::vector<double> &_firstValues) const { _firstValues = firstValues; };
    void getSecondValues(std::vector<double> &_secondValues) const { _secondValues = secondValues; };

    size_t getNumberWindowPoints() const { return numberWindowPoints; };
    size_t getNumberTestPoints() const { return numberTestPoints; };
    void getTestPoints(std::vector<point> &_testPoints) const { _testPoints = testPoints; };

    size_t getNumberCoefficients() const { return numberCoefficients; };
    void getCoefficients(std::vector<double> &_coefficients) const { _coefficients = coefficients; };

    void setGsaMethod(const GsaMethod<OneDimensionalSupportiveOptProblem> &_gsa) { gsa = _gsa; };

    void setIsSortX(bool _isSortX) { isSortX = _isSortX; };

    void getOptimalPoints(std::vector<std::vector<double>> &_optimalPoints) const override;

    double computeObjectiveFunction(const std::vector<double> &x) const override;
    double computeConstraintFunction(const std::vector<double> &x, size_t index) const override;
};

#endif // _BASE_FITTING_FAMILY_OPT_PROBLEMS_H_
