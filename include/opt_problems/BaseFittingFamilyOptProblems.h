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

    double alpha, delta;
    double leftBoundWideWindow, rightBoundWideWindow;
    std::vector<double> firstShortWindowValues, secondShortWindowValues;
    size_t numberWideWindowPoints, numberTestPoints;
    mutable std::vector<point> testPoints;

    size_t numberCoefficients;
    mutable std::vector<double> coefficients;

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
        double _alpha = 0.0, double _delta = 0.0,
        double _leftBoundWideWindow = 0.0, double _rightBoundWideWindow = 0.0,
        const std::vector<double> &insideWideWindowPoints = std::vector<double>{},
        double firstShortWindowPoint = 0.0, const std::vector<double> &_firstShortWindowValues = std::vector<double>{},
        double secondShortWindowPoint = 0.0, const std::vector<double> &_secondShortWindowValues = std::vector<double>{},
        double lastShortWindowPoint = 0.0,
        const GsaMethod<OneDimensionalSupportiveOptProblem> &_gsa = GsaMethod<OneDimensionalSupportiveOptProblem>(),
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
          gsa(_gsa),
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

    // TODO: implement any one-dimension opt method
    void setGsaMethod(const GsaMethod<OneDimensionalSupportiveOptProblem> &_gsa) { gsa = _gsa; };

    void setIsSortX(bool _isSortX) { isSortX = _isSortX; };

    void getOptimalPoints(std::vector<std::vector<double>> &_optimalPoints) const override;

    double computeObjectiveFunction(const std::vector<double> &x) const override;
    double computeConstraintFunction(const std::vector<double> &x, size_t index) const override;
};

#endif // _BASE_FITTING_FAMILY_OPT_PROBLEMS_H_
