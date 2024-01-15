#ifndef IMGO_METHOD_H_
#define IMGO_METHOD_H_

#include <vector>
#include <algorithm>
#include <limits>
#include <iterator>

#include <base_classes/opt_methods/IGeneralNumericalOptMethod.h>
#include <base_classes/opt_methods/ICharacteristicOptMethod.h>
#include <base_classes/opt_methods/IIndexSchemeOptMethod.h>
#include <base_structures/trials/IndexTrial.h>
#include <MyMath.h>

static const double epsilon = 1e-14;

template <typename OptProblemType>
class ImgoMethod : public opt::IGeneralNumericalOptMethod<opt::IndexTrial, OptProblemType>, 
    public opt::ICharacteristicOptMethod<opt::IndexTrial>, public opt::IIndexSchemeOptMethod
{
public:
    using GeneralNumericalMethod = opt::IGeneralNumericalOptMethod<opt::IndexTrial, OptProblemType>;
    using typename GeneralNumericalMethod::GeneralMethod;

    struct Parameters : public GeneralNumericalMethod::Parameters {
        std::vector<double> reliability;
        double d;

        Parameters(double _accuracy = 0.001, double _error = 0.001,
                   int _maxTrials = 1000, int _maxFevals = 1000,
                   std::vector<double> _reliability = std::vector<double>{}, double _d = 0.0):
            GeneralNumericalMethod::Parameters(_accuracy, _error, _maxTrials, _maxFevals),
            reliability(_reliability), d(_d)
        {};
    };

    using typename GeneralNumericalMethod::StoppingCondition;
    using typename GeneralNumericalMethod::StoppingConditions;

    struct Result : public GeneralNumericalMethod::Result {
        std::vector<double> constantsEstimation;

        Result(const typename OptProblemType::Point &_point = typename OptProblemType::Point(),
               double _value = 0.0, int _numberTrials = 0, int _numberFevals = 0,
               StoppingCondition _stoppingCondition = StoppingConditions::ACCURACY,
               std::vector<double> _constantsEstimation = std::vector<double>{})
            : GeneralNumericalMethod::Result(_point, _value, _numberTrials, _numberFevals, _stoppingCondition),
              constantsEstimation(_constantsEstimation) {};
    };

    class Report : public GeneralNumericalMethod::IReport {
    protected:
        void printPoint(std::ostream &stream, const typename OptProblemType::Point &point) const override;
        void printPoints(std::ostream &stream,
                         const std::vector<typename OptProblemType::Point> &points) const override;
        void printOptProblem(std::ostream &stream, const OptProblemType &optProblem) const override;
        void printMethodParameters(std::ostream &stream,
                                   const typename GeneralMethod::Parameters &parameters) const override;
        void printResultMethod(std::ostream &stream,
                               const typename GeneralMethod::Result &result) const override;
        void printErrorEstimate(std::ostream &stream, const OptProblemType &optProblem,
                                const typename GeneralMethod::Result &result) const override;

    public:
        Report():
            GeneralNumericalMethod::IReport()
        {};
    };

protected:
    double d;
    int lastI, lastTrialPosI;
    std::vector<bool> IIsCalc;

    void calcCharacteristic() override;
    void insertInSorted(const opt::IndexTrial &trial) override;

    opt::IndexTrial newTrial(const typename OptProblemType::Point &x) override;
    typename OptProblemType::Point selectNewPoint() override;
    double estimateSolution(typename OptProblemType::Point &x) const override;

    bool stopConditions() override;
    bool stopConditionsTest() override;

    void estimatingConstants() override;
    void calcZValues() override;

    void setResult(typename GeneralMethod::Result &result) const override;

public:
    ImgoMethod(const OptProblemType &_problem = OptProblemType(), const Parameters &parameters = Parameters())
        : GeneralNumericalMethod(_problem, parameters), opt::IIndexSchemeOptMethod(parameters.reliability),
          opt::ICharacteristicOptMethod<opt::IndexTrial>(), d(parameters.d), lastI(0), lastTrialPosI(0),
          IIsCalc(this->problem.getNumberConstraints() + 1) {}

    void setProblem(const OptProblemType &_problem) {
        GeneralMethod::setProblem(_problem);
        IIsCalc.resize(this->problem.getNumberConstraints() + 1);
    }

    void setParameters(const typename GeneralMethod::Parameters &parameters) override {
        GeneralNumericalMethod::setParameters(parameters);
        auto parametersCast = static_cast<const Parameters&>(parameters);
        setReliability(parametersCast.reliability);
        d = parametersCast.d;
    };
    void getParameters(typename GeneralMethod::Parameters &parameters) const override {
        parameters = Parameters(this->accuracy, this->error, this->maxTrials, this->maxFevals, reliability, d);
    };

    void setD(double _d) { d = _d; };
    double getD() const { return d; };

    typename GeneralNumericalMethod::Result* createResult() const override { return new Result(); };
    typename GeneralNumericalMethod::IReport* createReport() const override { return new Report(); };

    void solve(typename GeneralMethod::Result &result) override;
    bool solveTest(typename GeneralMethod::Result &result) override;
};

template <typename OptProblemType>
void ImgoMethod<OptProblemType>::Report::printPoint(
    std::ostream &stream, const typename OptProblemType::Point &point) const
{
    stream << "(" << point << ")";
}

template <typename OptProblemType>
void ImgoMethod<OptProblemType>::Report::printPoints(
    std::ostream &stream, const std::vector<typename OptProblemType::Point> &points) const
{
    stream << "[";
    printPoint(stream, points[0]);
    size_t numberPoints = points.size();
    for (size_t i = 1; i < numberPoints; ++i) {
        stream << "; ";
        printPoint(stream, points[i]);
    }
    stream << "]";
}

template <typename OptProblemType>
void ImgoMethod<OptProblemType>::Report::printOptProblem(
    std::ostream &stream, const OptProblemType &optProblem) const
{
    stream << "Number of constraints = " << optProblem.getNumberConstraints() << "\n";
    stream << "[a; b] = [" << optProblem.getSearchArea().lowerBound << "; " <<
                              optProblem.getSearchArea().upBound << "]"<< "\n";
    std::vector<double> optimalPoints;
    optProblem.getOptimalPoints(optimalPoints);
    if (!optimalPoints.empty()) {
        stream << "X* = ";
        printPoints(stream, optimalPoints);
        stream << "\n";
        stream << "f(X*) = " << optProblem.getOptimalValue() << "\n";
    }
    stream << "Lipschitz constants:" << "\n";
    std::vector<double> lipschitzConstants;
    optProblem.getLipschitzConstants(lipschitzConstants);
    size_t numberLipschitzConstants = lipschitzConstants.size();
    if (numberLipschitzConstants > 0) {
        stream << "L*(f) = " << lipschitzConstants[numberLipschitzConstants - 1] << "\n";
        for (size_t j = 0; j < numberLipschitzConstants - 1; ++j) {
            stream << "L*(g" << j + 1 << ") = " << lipschitzConstants[j] << "\n";
        }
    }
}

template <typename OptProblemType>
void ImgoMethod<OptProblemType>::Report::printMethodParameters(
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
}

template <typename OptProblemType>
void ImgoMethod<OptProblemType>::Report::printResultMethod(
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
void ImgoMethod<OptProblemType>::Report::printErrorEstimate(
    std::ostream &stream,
    const OptProblemType &optProblem,
    const typename GeneralMethod::Result &result) const
{
    double point = result.point;
    std::vector<double> optimalPoints;
    optProblem.getOptimalPoints(optimalPoints);
    if (!optimalPoints.empty()) {
        auto iter = std::min_element(optimalPoints.begin(), optimalPoints.end(),
            [&point] (const double &pointFirst, const double &pointSecond) {
                return std::abs(pointFirst - point) < std::abs(pointSecond - point);
            });
        stream << "|X* - X| = " << std::abs(*iter - point) << "\n";
        stream << "|f(X*) - f(X)| = " << std::abs(optProblem.getOptimalValue() - result.value) << "\n";
    }
}


template <typename OptProblemType>
void ImgoMethod<OptProblemType>::calcCharacteristic() {
    double R = -std::numeric_limits<double>::infinity(), Rtmp;
    double constantEstimation, zValue, dx;

    size_t trialPointsSize = this->trialPoints.size();
    for (size_t i = 1; i < trialPointsSize; ++i) {
        dx = this->trialPoints[i].x - this->trialPoints[i - 1].x;
        if (this->trialPoints[i].nu == this->trialPoints[i - 1].nu) {
            constantEstimation = constantsEstimation[this->trialPoints[i].nu];
            zValue = zValues[this->trialPoints[i].nu];
            Rtmp = dx + (this->trialPoints[i].z - this->trialPoints[i - 1].z) *
                        (this->trialPoints[i].z - this->trialPoints[i - 1].z) /
                        (reliability[this->trialPoints[i].nu] * reliability[this->trialPoints[i].nu] *
                         constantEstimation * constantEstimation * dx) -
                   2.0 * (this->trialPoints[i].z + this->trialPoints[i - 1].z - 2.0 * zValue) /
                   (reliability[this->trialPoints[i].nu]  * constantEstimation);
        } else if (this->trialPoints[i - 1].nu < this->trialPoints[i].nu) {
            constantEstimation = constantsEstimation[this->trialPoints[i].nu];
            zValue = zValues[this->trialPoints[i].nu];
            Rtmp = 2.0 * dx - 4.0 * (this->trialPoints[i].z - zValue) / 
                   (reliability[this->trialPoints[i].nu] * constantEstimation);
        } else {
            constantEstimation = constantsEstimation[this->trialPoints[i - 1].nu];
            zValue = zValues[this->trialPoints[i - 1].nu];
            Rtmp = 2.0 * dx - 4.0 * (this->trialPoints[i - 1].z - zValue) /
                   (reliability[this->trialPoints[i].nu] * constantEstimation);
        }
        if (Rtmp > R) {
            R = Rtmp;
            t = i;
        }
    }
}

template <typename OptProblemType>
void ImgoMethod<OptProblemType>::insertInSorted(const opt::IndexTrial &trial) {
    auto iter = this->trialPoints.begin();

    std::advance(iter, t);
    this->trialPoints.insert(iter, trial);

    iter = this->I[trial.nu].begin();
    auto iterEnd = this->I[trial.nu].end();

    while(true) {
        if (iter == iterEnd || iter->x > trial.x) break;
        iter++;
    }
    iter = this->I[trial.nu].insert(iter, trial);

    lastI = trial.nu;
    lastTrialPosI = distance(this->I[trial.nu].begin(), iter);

    // Binary search
    // iter = this->I[trial.nu - 1].begin();
    // iterEnd = this->I[trial.nu - 1].end();
    // int dist = this->I[trial.nu - 1].size();

    // if (dist != 0) {
    //     advance(iter, dist / 2);
    //     while (true) {
    //         dist /= 2;
    //         if (trial.x < iter->x) {
    //             advance(iter, -dist / 2);
    //         } else {
    //             if (iter  trial.x < (iter + 1)->x) {
    //                 break;
    //             } else {
    //                 advance(iter, dist / 2);
    //             }
    //         }
    //     }
    // }
    // iter = this->I[trial.nu - 1].insert(iter, trial);

    // lastI = trial.nu;
    // lastTrialPosI = distance(this->I[trial.nu - 1].begin(), iter);
}

template <typename OptProblemType>
void ImgoMethod<OptProblemType>::estimatingConstants() {
    // Constant complexity
    size_t sizeI = I[lastI].size();
    size_t numberFunctions = this->problem.getNumberConstraints() + 1;
    for (size_t nu = 0; nu < numberFunctions; ++nu) {
        if (!IIsCalc[nu]) constantsEstimation[nu] = 0.0;
    }
    if (sizeI >= 3) {
        if (lastTrialPosI == 0) {
            constantsEstimation[lastI] = std::max({ constantsEstimation[lastI],
                std::abs(I[lastI][1].z - I[lastI][0].z) / (I[lastI][1].x - I[lastI][0].x) });  
        } else if (lastTrialPosI == sizeI - 1) {
            constantsEstimation[lastI] = std::max({ constantsEstimation[lastI],
                std::abs(I[lastI][sizeI - 1].z - I[lastI][sizeI - 2].z) /
                        (I[lastI][sizeI - 1].x - I[lastI][sizeI - 2].x) });
        } else {
            constantsEstimation[lastI] = std::max({ constantsEstimation[lastI],
                std::abs(I[lastI][lastTrialPosI].z - I[lastI][lastTrialPosI - 1].z) / 
                        (I[lastI][lastTrialPosI].x - I[lastI][lastTrialPosI - 1].x),
                std::abs(I[lastI][lastTrialPosI + 1].z - I[lastI][lastTrialPosI].z) / 
                        (I[lastI][lastTrialPosI + 1].x - I[lastI][lastTrialPosI].x) });
        }
    } else if (sizeI == 2) {
        constantsEstimation[lastI] = std::max({ constantsEstimation[lastI],
            std::abs(I[lastI][1].z - I[lastI][0].z) / (I[lastI][1].x - I[lastI][0].x) });
    }
    if (std::abs(constantsEstimation[lastI]) > epsilon) IIsCalc[lastI] = true;
    for (size_t nu = 0; nu < numberFunctions; ++nu) {
        if (std::abs(constantsEstimation[nu]) <= epsilon) constantsEstimation[nu] = 1.0;
    }

    // with optimization(linear)
    // size_t nuLastTrial = lastTrial.nu;
    // size_t numberFunctions = this->problem.getNumberConstraints() + 1;
    // sizeI = I[nuLastTrial].size();
    // for (int nu = 0; nu < numberConstraints + 1; nu++) {
    //     if (!calcI[nuLastTrial]) mu[nu] = 0.0;
    // }
    // for (int i = 1; i < sizeI; i++) {
    //     muTmp = abs(I[nuLastTrial][i].z - I[nuLastTrial][(size_t)i - 1].z) /
    //                (I[nuLastTrial][i].x - I[nuLastTrial][(size_t)i - 1].x);
    //     if (muTmp > mu[nuLastTrial]) {
    //         mu[nuLastTrial] = muTmp;
    //         if (abs(mu[nuLastTrial]) > epsilon) calcI[nuLastTrial] = true;
    //     }
    // }
    // for (int nu = 0; nu < numberConstraints + 1; nu++) {
    //     if (abs(mu[nu]) <= epsilon) mu[nu] = 1.0;
    // }

    // without optimization
    // for (int nu = 0; nu < numberConstraints + 1; nu++) {
    //     mu[nu] = 0.0;
    // }
    // for (int nu = 0; nu < numberConstraints + 1; nu++) {
    //     sizeI = I[nu].size();
    //     for (int i = 1; i < sizeI; i++) {
    //         for (int j = 0; j < i; j++) {
    //             muTmp = abs(I[nu][i].z - I[nu][j].z) / (I[nu][i].x - I[nu][j].x);
    //             if (muTmp > mu[nu]) {
    //                 mu[nu] = muTmp;
    //             }
    //         }
    //     }
    //     if (abs(mu[nu]) <= epsilon) {
    //         mu[nu] = 1.0;
    //     };
    // }
}

template <typename OptProblemType>
void ImgoMethod<OptProblemType>::calcZValues() {
    size_t ISize;
    size_t numberFunctions = this->problem.getNumberConstraints() + 1;
    for (size_t nu = 0; nu < numberFunctions; ++nu) {
        if (I[nu].size() != 0) {
            zValues[nu] = I[nu][0].z;
            ISize = I[nu].size();
            for (size_t i = 1; i < ISize; ++i) {
                if (I[nu][i].z < zValues[nu]) {
                    zValues[nu] = I[nu][i].z;
                }
            }
            if (zValues[nu] <= 0.0 && nu != numberFunctions - 1) {
                zValues[nu] = -constantsEstimation[nu] * d;
            }
        }
    }
}

template <typename OptProblemType>
opt::IndexTrial ImgoMethod<OptProblemType>::newTrial(const typename OptProblemType::Point &x) {
    opt::IndexTrial trial(x);
    this->numberFevals++;
    double fValue;
    size_t numberConstraints = this->problem.getNumberConstraints();
    for (size_t j = 0; j < numberConstraints; ++j) {
        fValue = this->problem.computeConstraint(x, j);
        if (fValue > 0) {
            trial.z = fValue;
            trial.nu = j;
            return trial;
        }
    }
    trial.z = this->problem.computeConstraint(x, numberConstraints);
    trial.nu = numberConstraints;
    return trial;
}

template <typename OptProblemType>
typename OptProblemType::Point ImgoMethod<OptProblemType>::selectNewPoint() {
    estimatingConstants();
    calcZValues();
    calcCharacteristic();
    
    return this->trialPoints[t].nu != this->trialPoints[t - 1].nu ?
           (this->trialPoints[t].x + this->trialPoints[t - 1].x) / 2.0 :
           (this->trialPoints[t].x + this->trialPoints[t - 1].x) / 2.0 - 
           (this->trialPoints[t].z - this->trialPoints[t - 1].z) /
           (2.0 * reliability[this->trialPoints[t].nu] * constantsEstimation[this->trialPoints[t].nu]);
}

template <typename OptProblemType>
double ImgoMethod<OptProblemType>::estimateSolution(typename OptProblemType::Point &x) const {
    double z = std::numeric_limits<double>::infinity();
    x = 0.0;

    size_t sizeTrials = this->trialPoints.size();
    size_t numberConstraints = this->problem.getNumberConstraints();
    for (size_t i = 0; i < sizeTrials; ++i) {
        if (this->trialPoints[i].nu == numberConstraints && this->trialPoints[i].z < z) {
            z = this->trialPoints[i].z;
            x = this->trialPoints[i].x;
        }
    }

    return z;
}

template <typename OptProblemType>
bool ImgoMethod<OptProblemType>::stopConditions() {
    if (std::abs(this->trialPoints[t].x - this->trialPoints[t - 1].x) <= this->accuracy) {
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
bool ImgoMethod<OptProblemType>::stopConditionsTest() {
    std::vector<double> optimalPoints;
    this->problem.getOptimalPoints(optimalPoints);
    size_t numberOptimalPoints = optimalPoints.size();
    for (size_t i = 0; i < numberOptimalPoints; ++i) {
        if (std::abs(this->trialPoints[t].x - optimalPoints[i]) <= this->error) {
            this->stoppingCondition = StoppingConditions::ERROR;
            return true;
        }
    }
    return stopConditions();
}

template <typename OptProblemType>
void ImgoMethod<OptProblemType>::setResult(typename GeneralMethod::Result &result) const {
    GeneralNumericalMethod::setResult(result);

    auto& resultCast = static_cast<Result&>(result);
    resultCast.constantsEstimation = constantsEstimation;
}

template <typename OptProblemType>
void ImgoMethod<OptProblemType>::solve(typename GeneralMethod::Result &result) {
    size_t numberFunctions = this->problem.getNumberConstraints() + 1;
    for (size_t i = 0; i < numberFunctions; ++i) {
        this->I[i].clear();
        IIsCalc[i] = false;
    }
    this->trialPoints.clear();
    this->numberFevals = 0;

    t = 0;
    insertInSorted(newTrial(this->problem.getSearchArea().lowerBound));
    t = 1;
    insertInSorted(newTrial(this->problem.getSearchArea().upBound));
    
    this->numberTrials = 2;

    opt::IndexTrial trial;
    double xNew;
    while(!stopConditions()) {
        xNew = selectNewPoint();

        trial = newTrial(xNew);
        this->numberTrials++;

        insertInSorted(trial);
    }

    setResult(result);
}

template <typename OptProblemType>
bool ImgoMethod<OptProblemType>::solveTest(typename GeneralMethod::Result &result) {
    size_t numberFunctions = this->problem.getNumberConstraints() + 1;
    for (size_t i = 0; i < numberFunctions; ++i) {
        this->I[i].clear();
        IIsCalc[i] = false;
    }
    this->trialPoints.clear();
    this->numberFevals = 0;

    insertInSorted(newTrial(this->problem.getSearchArea().lowerBound));
    insertInSorted(newTrial(this->problem.getSearchArea().upBound));
    t = 1;
    
    this->numberTrials = 2;

    opt::IndexTrial trial;
    double xNew;
    while(!stopConditionsTest()) {
        xNew = selectNewPoint();

        trial = newTrial(xNew);
        this->numberTrials++;

        insertInSorted(trial);
    }

    setResult(result);

    return static_cast<Result&>(result).stoppingCondition == StoppingConditions::ERROR;
}

#endif // IMGO_METHOD_H_
