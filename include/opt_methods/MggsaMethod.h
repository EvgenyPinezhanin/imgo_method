#ifndef _MGGSA_METHOD_H
#define _MGGSA_METHOD_H

#include <vector>
#include <algorithm>
#include <limits>
#include <iterator>

#include <general/classes/opt_methods/IGeneralNumericalOptMethod.h>
#include <general/classes/opt_methods/ICharacteristicOptMethod.h>
#include <general/classes/opt_methods/IIndexSchemeOptMethod.h>
#include <general/structures/trials/IndexTrial.h>
#include <general/structures/search_areas/MultiDimensionalSearchArea.h>
#include <MyMath.h>
#include <Map.h>

static const double peanoA = 0.0, peanoB = 1.0, peanoRandom = 0.5;

template <typename OptProblemType>
class MggsaMethod : public opt::IGeneralNumericalOptMethod<opt::IndexTrial, OptProblemType>, 
    public opt::ICharacteristicOptMethod<opt::IndexTrial>, public opt::IIndexSchemeOptMethod
{
public:
    using GeneralNumericalMethod = opt::IGeneralNumericalOptMethod<opt::IndexTrial, OptProblemType>;
    using typename GeneralNumericalMethod::GeneralMethod;

    enum class TypeSolve { SOLVE, RESOLVE };

    struct Parameters : public GeneralNumericalMethod::Parameters {
        std::vector<double> reliability;
        double d;
        size_t density, key, increment;
        TypeSolve typeSolve;

        Parameters(double _accuracy = 0.001, double _error = 0.001,
                   int _maxTrials = 1000, int _maxFevals = 1000,
                   std::vector<double> _reliability = std::vector<double>{}, double _d = 0.0,
                   size_t _density = 10, size_t _key = 1, size_t _increment = 0,
                   TypeSolve _typeSolve = TypeSolve::SOLVE)
            : GeneralNumericalMethod::Parameters(_accuracy, _error, _maxTrials, _maxFevals),
              reliability(_reliability), d(_d), density(_density), key(_key),
              increment(_increment), typeSolve(_typeSolve) {};
    };

    using typename GeneralNumericalMethod::StoppingCondition;

    struct StoppingConditions : public GeneralNumericalMethod::StoppingConditions {
        static const StoppingCondition COINCIDENCE = 4;
    };

    struct Result : public GeneralNumericalMethod::Result {
        std::vector<double> constantsEstimation;

        Result(const typename OptProblemType::Point &_point = typename OptProblemType::Point(),
               double _value = 0.0, int _numberTrials = 0, int _numberFevals = 0,
               double _resultingAccuracy = 0.0,
               StoppingCondition _stoppingCondition = StoppingConditions::ACCURACY,
               std::vector<double> _constantsEstimation = std::vector<double>{})
            : GeneralNumericalMethod::Result(_point, _value, _numberTrials, _numberFevals,
              _resultingAccuracy, _stoppingCondition), constantsEstimation(_constantsEstimation) {};
    };

    class Report : public GeneralNumericalMethod::IReport {
    protected:
        void printOptProblem(std::ostream &stream, const OptProblemType &optProblem) const override;
        void printMethodParameters(std::ostream &stream,
                                   const typename GeneralMethod::Parameters &parameters) const override;
        void printResultMethod(std::ostream &stream,
                               const typename GeneralMethod::Result &result) const override;
        void printErrorEstimate(std::ostream &stream, const OptProblemType &optProblem,
                                const typename GeneralMethod::Result &result) const override;

    public:
        Report() : GeneralNumericalMethod::IReport() {};

        void printPoint(std::ostream &stream, const typename OptProblemType::Point &point) const override;
    };

protected:
    double d;
    size_t density, key, increment;

    size_t lastI, lastTrialPosI;
    std::vector<size_t> lastTrialsPosI;
    std::vector<bool> IIsCalc;
    TypeSolve typeSolve;

    std::vector<std::vector<double>> points;

    double newOneDimensionX, h;
    size_t M;

    void calcCharacteristic() override;
    void insertInSorted(const opt::IndexTrial &trial) override;

    opt::IndexTrial newTrial(const typename OptProblemType::Point &x) override;
    typename OptProblemType::Point selectNewPoint() override;
    double estimateSolution(typename OptProblemType::Point &x) const override;

    void estimatingConstants() override;
    void calcZValues() override;

    double calcH(double xNew);
    bool checkDensity(double h);
    void x(const std::vector<double> &P, std::vector<double> &X);

    bool stopConditions() override;
    bool stopConditionsTest() override;

    void setResult(typename GeneralMethod::Result &result) const override;

public:
    MggsaMethod(const OptProblemType &_problem = OptProblemType(), const Parameters &parameters = Parameters())
        : GeneralNumericalMethod(_problem, parameters), opt::IIndexSchemeOptMethod(parameters.reliability),
          opt::ICharacteristicOptMethod<opt::IndexTrial>(), d(parameters.d), density(parameters.density),
          key(parameters.key), increment(parameters.increment), typeSolve(parameters.typeSolve), lastI(),
          lastTrialsPosI(), IIsCalc(this->problem.getNumberConstraints() + 1), points(), newOneDimensionX(0.0),
          h(0.0), M(0) {}

    void setProblem(const OptProblemType &_problem) {
        GeneralMethod::setProblem(_problem);
        IIsCalc.resize(this->problem.getNumberConstraints() + 1);
    }

    void setParameters(const typename GeneralMethod::Parameters &parameters) override {
        GeneralNumericalMethod::setParameters(parameters);
        auto parametersCast = static_cast<const Parameters&>(parameters);
        setReliability(parametersCast.reliability);
        d = parametersCast.d;
        density = parametersCast.density;
        key = parametersCast.key;
        increment = parametersCast.increment;
        typeSolve = parametersCast.typeSolve;
    };
    void getParameters(typename GeneralMethod::Parameters &parameters) const override {
        parameters = Parameters(this->accuracy, this->error, this->maxTrials, this->maxFevals, reliability,
                                d, density, key, increment, typeSolve);
    };

    void y(double x, std::vector<double> &X) const;

    void setD(double _d) { d = _d; };
    double getD() const { return d; };

    void setDensity(size_t _density) { density = _density; };
    size_t getDensity() const { return density; };

    void setKey(size_t _key) { key = _key; };
    size_t getKey() const { return key; };

    void setIncrement(size_t _increment) { increment = _increment; };
    size_t getIncrement() { return increment; };

    void setTypeSolve(TypeSolve _typeSolve) { typeSolve = _typeSolve; };
    TypeSolve getTypeSolve() const { return typeSolve; };

    void getPoints(std::vector<std::vector<double>> &points) const { points = this->points; };

    typename GeneralNumericalMethod::Result* createResult() const override { return new Result(); };
    typename GeneralNumericalMethod::IReport* createReport() const override { return new Report(); };

    void solve(typename GeneralMethod::Result &result) override;
    bool solveTest(typename GeneralMethod::Result &result) override;
};

template <typename OptProblemType>
void MggsaMethod<OptProblemType>::Report::printOptProblem(
    std::ostream &stream, const OptProblemType &optProblem) const
{
    stream << "Number of constraints = " << optProblem.getNumberConstraints() << "\n";
    stream << "[a; b] = [";
    printPoint(stream, optProblem.getSearchArea().lowerBound);
    stream << "; ";
    printPoint(stream, optProblem.getSearchArea().upBound);
    stream << "]"<< "\n";
    std::vector<std::vector<double>> optimalPoints;
    optProblem.getOptimalPoints(optimalPoints);
    if (!optimalPoints.empty()) {
        stream << "X* = ";
        this->printPoints(stream, optimalPoints);
        stream << "\n";
        stream << "f(X*) = " << optProblem.getOptimalValue() << "\n";
    }
    double objectiveLipschitzConstant = optProblem.getObjectiveLipschitzConstant();
    if (objectiveLipschitzConstant != 0.0) {
        stream << "Lipschitz constants:" << "\n";
        stream << "L*(f) = " << objectiveLipschitzConstant << "\n";
        std::vector<double> constraintLipschitzConstants;
        optProblem.getConstraintLipschitzConstants(constraintLipschitzConstants);
        if (!constraintLipschitzConstants.empty()) {
            size_t numberLipschitzConstants = optProblem.getNumberConstraints();
            for (size_t j = 0; j < numberLipschitzConstants; ++j) {
                stream << "L*(g" << j + 1 << ") = " << constraintLipschitzConstants[j] << "\n";
            }
        }
    }
}

template <typename OptProblemType>
void MggsaMethod<OptProblemType>::Report::printMethodParameters(
    std::ostream &stream, const typename GeneralMethod::Parameters &parameters) const
{
    GeneralNumericalMethod::IReport::printMethodParameters(stream, parameters);
    auto parametersCast = static_cast<const Parameters&>(parameters);

    stream << "reliability = (";
    size_t numberConstraints = parametersCast.reliability.size();
    for (size_t i = 0; i < numberConstraints - 1; ++i) {
        stream << parametersCast.reliability[i] << ", ";
    }
    stream << parametersCast.reliability[numberConstraints - 1] << ")" << "\n";
    stream << "d = " << parametersCast.d << "\n";
    stream << "density = " << parametersCast.density << "\n";
    stream << "key = " << parametersCast.key << "\n";
    stream << "increment = " << parametersCast.increment << "\n";
    // TODO: add typeSolve
}

template <typename OptProblemType>
void MggsaMethod<OptProblemType>::Report::printResultMethod(
    std::ostream &stream,
    const typename GeneralMethod::Result &result) const
{
    GeneralNumericalMethod::IReport::printResultMethod(stream, result);
    
    auto resultCast = static_cast<const Result&>(result);
    stream << "Estimation of the Lipschitz constant:" << "\n";
    size_t numberLipschitzConstants = resultCast.constantsEstimation.size();
    stream << "L(f) = " << resultCast.constantsEstimation[numberLipschitzConstants - 1] << "\n";
    for (size_t j = 0; j < numberLipschitzConstants - 1; ++j) {
        stream << "L(g" << j + 1 << ") = " << resultCast.constantsEstimation[j] << "\n";
    }
}

template <typename OptProblemType>
void MggsaMethod<OptProblemType>::Report::printErrorEstimate(
    std::ostream &stream,
    const OptProblemType &optProblem,
    const typename GeneralMethod::Result &result) const
{
    std::vector<double> point = result.point;
    std::vector<std::vector<double>> optimalPoints;
    optProblem.getOptimalPoints(optimalPoints);
    if (!optimalPoints.empty()) {
        auto iter = std::min_element(optimalPoints.begin(), optimalPoints.end(),
            [&point] (const std::vector<double> &pointFirst, const std::vector<double> &pointSecond) {
                return euclideanDistance(pointFirst, point) < euclideanDistance(pointSecond, point);
            });
        stream << "||X* - X|| = " << euclideanDistance(*iter, point) << "\n";
        stream << "|f(X*) - f(X)| = " << std::abs(optProblem.getOptimalValue() - result.value) << "\n";
    }
}

template <typename OptProblemType>
void MggsaMethod<OptProblemType>::Report::printPoint(
    std::ostream &stream, const typename OptProblemType::Point &point) const
{
    stream << "(" << point[0];
    size_t dimensionalPoint = point.size();
    for (size_t i = 1; i < dimensionalPoint; ++i) {
        stream << ", " << point[i];
    }
    stream << ")";
}


template <typename OptProblemType>
void MggsaMethod<OptProblemType>::calcCharacteristic() {
    double R = -std::numeric_limits<double>::infinity(), Rtmp;
    double constantEstimation, zValue, dx;
    opt::IndexTrial lastTrial = this->trialPoints[0], nextTrial;
    size_t dimension = this->problem.getSearchArea().dimension;

    size_t sizeTrialPoints = this->trialPoints.size();
    for (size_t i = 1; i < sizeTrialPoints; i++) {
        nextTrial = this->trialPoints[i];
        dx = std::pow(nextTrial.x - lastTrial.x, 1.0 / dimension);
        if (nextTrial.nu == lastTrial.nu) {
            constantEstimation = constantsEstimation[nextTrial.nu];
            zValue = zValues[nextTrial.nu];
            Rtmp = dx + (nextTrial.z - lastTrial.z) * (nextTrial.z - lastTrial.z) /
                   (reliability[this->trialPoints[i].nu] * reliability[this->trialPoints[i].nu]*
                   constantEstimation * constantEstimation * dx) - 2.0 * (nextTrial.z + lastTrial.z - 2.0 * zValue) /
                   (reliability[this->trialPoints[i].nu] * constantEstimation);
        } else if (lastTrial.nu < nextTrial.nu) {
            constantEstimation = constantsEstimation[nextTrial.nu];
            zValue = zValues[nextTrial.nu];
            Rtmp = 2.0 * dx - 4.0 * (nextTrial.z - zValue) /
                   (reliability[nextTrial.nu] * constantEstimation);
        } else  {
            constantEstimation = constantsEstimation[lastTrial.nu];
            zValue = zValues[lastTrial.nu];
            Rtmp = 2.0 * dx - 4.0 * (lastTrial.z - zValue) /
                   (reliability[lastTrial.nu] * constantEstimation);
        }
        if (Rtmp > R) {
            R = Rtmp;
            t = i;
        }
        lastTrial = nextTrial;
    }
}

template <typename OptProblemType>
void MggsaMethod<OptProblemType>::insertInSorted(const opt::IndexTrial &trial) {
    auto iter = this->trialPoints.begin();

    if (key == 3) {
        auto iterEnd = this->trialPoints.end();

        while(true) {
            if (iter == iterEnd || iter->x > trial.x) break;
            iter++;
        }
        iter = this->trialPoints.insert(iter, trial);
    } else {
        std::advance(iter, t);
        this->trialPoints.insert(iter, trial);
    }

    iter = this->I[trial.nu].begin();
    auto iterEnd = this->I[trial.nu].end();

    while(true) {
        if (iter == iterEnd || iter->x > trial.x) break;
        iter++;
    }
    iter = this->I[trial.nu].insert(iter, trial);

    lastI = trial.nu;
    lastTrialPosI = distance(this->I[trial.nu].begin(), iter);
    // TODO: binary search
}

template <typename OptProblemType>
void MggsaMethod<OptProblemType>::estimatingConstants() {
    // Linear complexity
    double constantEstimationTmp;
    size_t sizeLastTrials = lastTrialsPosI.size();
    size_t sizeI = I[lastI].size();
    size_t dimension = this->problem.getSearchArea().dimension;

    size_t numberFunctions = this->problem.getNumberConstraints() + 1;
    for (size_t nu = 0; nu < numberFunctions; ++nu) {
        if (!IIsCalc[nu]) constantsEstimation[nu] = 0.0;
    }
    for (size_t i = 0; i < sizeLastTrials; ++i) {
        for (size_t j = 0; j < sizeI; ++j) {
            if (I[lastI][j].x != I[lastI][lastTrialsPosI[i]].x) {
                constantEstimationTmp = std::abs(I[lastI][j].z - I[lastI][lastTrialsPosI[i]].z) /
                                        std::pow(std::abs(I[lastI][j].x - I[lastI][lastTrialsPosI[i]].x), 1.0 / dimension);
                if (constantEstimationTmp > constantsEstimation[lastI]) {
                    constantsEstimation[lastI] = constantEstimationTmp;
                    if (std::abs(constantsEstimation[lastI]) > this->epsilon) IIsCalc[lastI] = true;
                }
            }
        }
    }
    for (size_t nu = 0; nu < numberFunctions; ++nu) {
        if (std::abs(constantsEstimation[nu]) <= this->epsilon) constantsEstimation[nu] = 1.0;
    }

    // TODO: recheck this code
    // with optimization (const) (not work)
    // nuLastTrials = lastTrials[0].nu - 1;
    // sizeI = I[nuLastTrials].size();
    // sizeLastTrials = lastTrials.size();
    // for (int nu = 0; nu < numberConstraints + 1; nu++) {
    //     if (!calcI[nu]) mu[nu] = 0.0;
    // }
    // if (I[nuLastTrials].size() >= 3) {
    //     for (int i = 0; i < sizeLastTrials; i++) {
    //         if (lastTrialsPos[i] == 0) {
    //             mu[nuLastTrials] = max({ mu[nuLastTrials], abs(I[nuLastTrials][1].z - I[nuLastTrials][0].z) /
    //                                                        pow(I[nuLastTrials][1].x - I[nuLastTrials][0].x, 1.0 / n) });  
    //         } else if (lastTrialsPos[i] == sizeI - 1) {
    //             mu[nuLastTrials] = max({ mu[nuLastTrials],
    //                                      abs(I[nuLastTrials][sizeI - 1].z - I[nuLastTrials][sizeI - 2].z) / 
    //                                      pow(I[nuLastTrials][sizeI - 1].x - I[nuLastTrials][sizeI - 2].x, 1.0 / n) });
    //         } else {
    //             mu[nuLastTrials] = max({ mu[nuLastTrials],
    //                                      abs(I[nuLastTrials][lastTrialsPos[i]].z - I[nuLastTrials][lastTrialsPos[i] - 1].z) / 
    //                                      pow(I[nuLastTrials][lastTrialsPos[i]].x - I[nuLastTrials][lastTrialsPos[i] - 1].x, 1.0 / n),
    //                                      abs(I[nuLastTrials][lastTrialsPos[i] + 1].z - I[nuLastTrials][lastTrialsPos[i]].z) / 
    //                                      pow(I[nuLastTrials][lastTrialsPos[i] + 1].x - I[nuLastTrials][lastTrialsPos[i]].x, 1.0 / n) });
    //         }
    //     }
    // } else if (I[nuLastTrials].size() == 2) {
    //     mu[nuLastTrials] = max({ mu[nuLastTrials], abs(I[nuLastTrials][1].z - I[nuLastTrials][0].z) /
    //                                                pow(I[nuLastTrials][1].x - I[nuLastTrials][0].x, 1.0 / n) });
    // }
    // if (mu[nuLastTrials] > epsilon) calcI[nuLastTrials] = true;
    // for (int nu = 0; nu < numberConstraints + 1; nu++) {
    //     if (mu[nu] <= epsilon) mu[nu] = 1.0;
    // }

    // without optimization
    // for (int nu = 0; nu < numberConstraints + 1; nu++) {
    //     mu[nu] = 0.0;
    // }
    // for (int nu = 0; nu < numberConstraints + 1; nu++) {
    //     sizeI = I[nu].size();
    //     for (int i = 1; i < sizeI; i++) { // при i = 0 - нет j
    //         for (int j = 0; j < i; j++) {
    //             muTmp = std::abs(I[nu][i].z - I[nu][j].z) / pow(I[nu][i].x - I[nu][j].x, 1.0 / dimension);
    //             if (muTmp > mu[nu]) {
    //                 mu[nu] = muTmp;
    //             }
    //         }
    //     }
    //     if (std::abs(mu[nu]) <= epsilon) {
    //         mu[nu] = 1.0;
    //     };
    // }
}

template <typename OptProblemType>
void MggsaMethod<OptProblemType>::calcZValues() {
    size_t sizeI, numberFunctions = this->problem.getNumberConstraints() + 1;
    for (size_t nu = 0; nu < numberFunctions; ++nu) {
        if (I[nu].size() != 0) {
            if (nu != M) {
                zValues[nu] = -constantsEstimation[nu] * d;
            } else {
                zValues[nu] = I[nu][0].z;
                sizeI = I[nu].size();
                for (size_t i = 1; i < sizeI; ++i) {
                    if (I[nu][i].z < zValues[nu]) {
                        zValues[nu] = I[nu][i].z;
                    }
                }
            }
        }
    }
}

template <typename OptProblemType>
opt::IndexTrial MggsaMethod<OptProblemType>::newTrial(const typename OptProblemType::Point &x) {
    opt::IndexTrial trial(h);
    ++this->numberFevals;
    double fValue;
    size_t numberConstraints = this->problem.getNumberConstraints();
    for (size_t j = 0; j < numberConstraints; ++j) {
        fValue = this->problem.computeConstraintFunction(x, j);
        if (fValue > 0) {
            trial.z = fValue;
            trial.nu = j;
            return trial;
        }
    }
    trial.z = this->problem.computeObjectiveFunction(x);
    trial.nu = numberConstraints;
    return trial;
}

template <typename OptProblemType>
typename OptProblemType::Point MggsaMethod<OptProblemType>::selectNewPoint() {
    estimatingConstants();
    calcZValues();
    calcCharacteristic();

    size_t dimension = this->problem.getSearchArea().dimension;
    newOneDimensionX = this->trialPoints[t].nu != this->trialPoints[t - 1].nu ?
                         (this->trialPoints[t].x + this->trialPoints[t - 1].x) / 2.0 :
                         (this->trialPoints[t].x + this->trialPoints[t - 1].x) / 2.0 -
                         sgn(this->trialPoints[t].z - this->trialPoints[t - 1].z) /
                         (2.0 * constantsEstimation[this->trialPoints[t].nu]) *
                         std::pow(std::abs(this->trialPoints[t].z - this->trialPoints[t - 1].z) /
                         constantsEstimation[this->trialPoints[t].nu], dimension);

    typename OptProblemType::Point X(dimension);
    y(newOneDimensionX, X);

    return X;
}

template <typename OptProblemType>
double MggsaMethod<OptProblemType>::estimateSolution(typename OptProblemType::Point &x) const {
    double z = std::numeric_limits<double>::infinity(), oneDimensionX = 0.0;

    size_t sizeTrials = this->trialPoints.size();
    size_t numberConstraints = this->problem.getNumberConstraints();
    for (size_t i = 0; i < sizeTrials; ++i) {
        if (this->trialPoints[i].nu == numberConstraints && this->trialPoints[i].z < z) {
            z = this->trialPoints[i].z;
            oneDimensionX = this->trialPoints[i].x;
        }
    }
    y(oneDimensionX, x);

    return z;
}

template <typename OptProblemType>
double MggsaMethod<OptProblemType>::calcH(double xNew) {
    size_t dimension = this->problem.getSearchArea().dimension;
    double d = 1.0 / (std::pow(2.0, density * dimension) * (std::pow(2.0, dimension - 1.0)));
    double h = floor(xNew / d) * d;
    return h;
}

template <typename OptProblemType>
bool MggsaMethod<OptProblemType>::checkDensity(double h) {
    size_t sizeTrialPoints = this->trialPoints.size();
    for (int i = 0; i < sizeTrialPoints; i++) {
        if (std::abs(h - this->trialPoints[i].x) <= this->epsilon) return false; 
    }
    return true;
}

template <typename OptProblemType>
void MggsaMethod<OptProblemType>::y(double x, std::vector<double> &X) const {
    int d = (key != 3) ? density : density + 1;
    size_t dimension = this->problem.getSearchArea().dimension;
    X.resize(dimension);
    mapd(x, d, X.data(), dimension, key);

    opt::MultiDimensionalSearchArea searchArea = this->problem.getSearchArea();
    for (int i = 0; i < dimension; i++) {
        X[i] = X[i] * (searchArea.upBound[i] - searchArea.lowerBound[i]) +
               (searchArea.lowerBound[i] + searchArea.upBound[i]) / 2.0;
    }
}

template <typename OptProblemType>
void MggsaMethod<OptProblemType>::x(const std::vector<double> &P, std::vector<double> &X) {
    // TODO: add dimension in class
    size_t dimension = this->problem.getSearchArea().dimension;
    std::vector<double> pCorrect(dimension);
    int sizeX;
    X.resize((size_t)pow(2, dimension));
    opt::MultiDimensionalSearchArea searchArea = this->problem.getSearchArea();
    for (int i = 0; i < dimension; i++) {
        pCorrect[i] = (P[i] - (searchArea.lowerBound[i] + searchArea.upBound[i]) / 2.0) /
                      (searchArea.upBound[i] - searchArea.lowerBound[i]);
    }

    invmad(density + 1, X.data(), (int)pow(2, dimension), &sizeX, pCorrect.data(), dimension, increment);
    X.resize(sizeX);
}

template <typename OptProblemType>
bool MggsaMethod<OptProblemType>::stopConditions() {
    size_t dimension = this->problem.getSearchArea().dimension;
    if (std::pow(this->trialPoints[t].x - this->trialPoints[t - 1].x, 1.0 / dimension) <= this->accuracy) {
        this->stoppingCondition = StoppingConditions::ACCURACY;
        return true;
    }
    if (this->numberTrials >= this->maxTrials) {
        this->stoppingCondition = StoppingConditions::MAXTRIALS;
        return true;
    }
    if (this->numberFevals >= this->maxFevals) {
        this->stoppingCondition = StoppingConditions::MAXFEVALS;
        return true;
    }
    return false;
}

template <typename OptProblemType>
bool MggsaMethod<OptProblemType>::stopConditionsTest() {
    std::vector<std::vector<double>> optimalPoints;
    this->problem.getOptimalPoints(optimalPoints);
    std::vector<double> X;
    y(this->trialPoints[t].x, X);

    size_t numberOptimalPoints = optimalPoints.size();
    for (size_t i = 0; i < numberOptimalPoints; ++i) {
        if (euclideanDistance(X, optimalPoints[i]) <= this->error) {
            this->stoppingCondition = StoppingConditions::ERROR;
            return true;
        }
    }
    return stopConditions();
}

template <typename OptProblemType>
void MggsaMethod<OptProblemType>::setResult(typename GeneralMethod::Result &result) const {
    GeneralNumericalMethod::setResult(result);

    auto& resultCast = static_cast<Result&>(result);
    resultCast.constantsEstimation = constantsEstimation;
}

template <typename OptProblemType>
void MggsaMethod<OptProblemType>::solve(typename GeneralMethod::Result &result) {
    opt::IndexTrial trial;
    std::vector<double> hNu, P, point;

    if (typeSolve == TypeSolve::SOLVE) {
        size_t numberConstraints = this->problem.getNumberConstraints();
        for (size_t nu = 0; nu < numberConstraints + 1; ++nu) {
            I[nu].clear();
            IIsCalc[nu] = false;
        }
        points.clear();
        this->trialPoints.clear();
        this->numberFevals = 0;

        this->trialPoints.push_back(opt::IndexTrial{ peanoA, -1.0, -1 });
        this->trialPoints.push_back(opt::IndexTrial{ peanoB, -1.0, -1 });
        t = 1;
        lastI = 0;
        lastTrialsPosI = std::vector<size_t>{ 1 };

        lastTrialsPosI.clear();
        if (key == 3) {
            h = calcH(peanoRandom);
            y(h, P);
            x(P, hNu);

            trial = newTrial(P);
            size_t sizeHNu = hNu.size();
            for (size_t i = 0; i < sizeHNu; ++i) {
                trial.x = hNu[i];
                insertInSorted(trial);
                lastTrialsPosI.push_back(lastTrialPosI);
            }
        } else {
            y(peanoRandom, P);
            h = peanoRandom;

            trial = newTrial(P);
            insertInSorted(trial);
            lastTrialsPosI.push_back(lastTrialPosI);
        }

        point = P;
        point.push_back(trial.z);
        points.push_back(point);

        this->numberTrials = 1;
        M = lastI;
    }

    while(!stopConditions()) {
        selectNewPoint();

        lastTrialsPosI.clear();
        if (key == 3) {
            h = calcH(newOneDimensionX);
            y(h, P);
            x(P, hNu);

            if (!checkDensity(h)) {
                this->stoppingCondition = StoppingConditions::COINCIDENCE;
                break;
            }

            trial = newTrial(P);
            size_t sizeHNu = hNu.size();
            for (size_t i = 0; i < sizeHNu; ++i) {
                trial.x = hNu[i];
                insertInSorted(trial);
                lastTrialsPosI.push_back(lastTrialPosI);
            }
        } else {
            h = newOneDimensionX;
            y(newOneDimensionX, P);
            trial = newTrial(P);
            insertInSorted(trial);
            lastTrialsPosI.push_back(lastTrialPosI);

        }
        ++this->numberTrials;

        point = P;
        point.push_back(trial.z);
        points.push_back(point);

        if (lastI > M) {
            M = lastI;
        }
    }
    
    setResult(result);
}

template <typename OptProblemType>
bool MggsaMethod<OptProblemType>::solveTest(typename GeneralMethod::Result &result) {
    opt::IndexTrial trial;
    std::vector<double> hNu, P, point;

    if (typeSolve == TypeSolve::SOLVE) {
        size_t numberConstraints = this->problem.getNumberConstraints();
        for (size_t nu = 0; nu < numberConstraints + 1; ++nu) {
            I[nu].clear();
            IIsCalc[nu] = false;
        }
        points.clear();
        this->trialPoints.clear();
        this->numberFevals = 0;

        this->trialPoints.push_back(opt::IndexTrial{ peanoA, -1.0, -1 });
        this->trialPoints.push_back(opt::IndexTrial{ peanoB, -1.0, -1 });
        t = 1;
        lastI = 0;
        lastTrialsPosI = std::vector<size_t>{ 1 };

        lastTrialsPosI.clear();
        if (key == 3) {
            h = calcH(peanoRandom);
            y(h, P);
            x(P, hNu);

            trial = newTrial(P);
            size_t sizeHNu = hNu.size();
            for (size_t i = 0; i < sizeHNu; ++i) {
                trial.x = hNu[i];
                insertInSorted(trial);
                lastTrialsPosI.push_back(lastTrialPosI);
            }
        } else {
            y(peanoRandom, P);
            h = peanoRandom;

            trial = newTrial(P);
            insertInSorted(trial);
            lastTrialsPosI.push_back(lastTrialPosI);
        }

        point = P;
        point.push_back(trial.z);
        points.push_back(point);

        this->numberTrials = 1;
        M = lastI;
    }

    while(!stopConditionsTest()) {
        selectNewPoint();

        lastTrialsPosI.clear();
        if (key == 3) {
            h = calcH(newOneDimensionX);
            y(h, P);
            x(P, hNu);

            if (!checkDensity(h)) {
                this->stoppingCondition = StoppingConditions::COINCIDENCE;
                break;
            }

            trial = newTrial(P);
            size_t sizeHNu = hNu.size();
            for (size_t i = 0; i < sizeHNu; ++i) {
                trial.x = hNu[i];
                insertInSorted(trial);
                lastTrialsPosI.push_back(lastTrialPosI);
            }
        } else {
            h = newOneDimensionX;
            y(newOneDimensionX, P);
            trial = newTrial(P);
            insertInSorted(trial);
            lastTrialsPosI.push_back(lastTrialPosI);

        }
        ++this->numberTrials;

        point = P;
        point.push_back(trial.z);
        points.push_back(point);

        if (lastI > M) {
            M = lastI;
        }
    }

    setResult(result);

    return static_cast<Result&>(result).stoppingCondition == StoppingConditions::ERROR;
}

#endif // _MGGSA_METHOD_H
