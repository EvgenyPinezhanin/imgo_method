#ifndef SCANNING_METHOD_H_
#define SCANNING_METHOD_H_

#include <algorithm>
#include <limits>
#include <iterator>

#include <base_classes/opt_methods/IGeneralNumericalOptMethod.h>
#include <base_classes/opt_methods/ICharacteristicOptMethod.h>
#include <opt_problems/OneDimensionalProblem.h>
#include <trials/Trial.h>
#include <my_math.h>

template <typename OptProblemType>
class ScanningMethod : public opt::ICharacteristicOptMethod<Trial>,
    public opt::IGeneralNumericalOptMethod<Trial, OptProblemType>
{
public:
    using GeneralNumMethod = opt::IGeneralNumericalOptMethod<Trial, OptProblemType>;
    using typename GeneralNumMethod::GeneralMethod;

    using typename GeneralNumMethod::Parameters;
    using typename GeneralNumMethod::StoppingConditions;
    using typename GeneralNumMethod::Result;

    struct Task : public GeneralNumMethod::Task {
        std::string blockName;
        int functionNumber;

        Task(const std::string &_name, std::string _blockName, int _functionNumber,
             const OptProblemType &_problem, Parameters &_parameters, bool _use = true):
            GeneralNumMethod::Task(_name, _problem, _parameters, _use),
            blockName(_blockName),
            functionNumber(_functionNumber)
        {};
    };

    class Report : public GeneralNumMethod::IReport {
    protected:
        void printOptProblem(std::ostream &stream, const OptProblemType &optProblem) const override;
        void printErrorEstimate(std::ostream &stream, const OptProblemType &optProblem,
                                const typename GeneralMethod::Result &result) const override;

    public:
        Report():
            GeneralNumMethod::IReport()
        {};
    };
private:
    int cycleBoundary;

protected:
    void insertInSorted(const Trial &trial) override;

    void calcCharacteristic() override;

    Trial newTrial(const typename OptProblemType::Point &x) override;
    typename OptProblemType::Point selectNewPoint() override;

    double estimateSolution(typename OptProblemType::Point &x) const override;

    bool stopConditions() override;
    bool stopConditionsTest() override;

    using GeneralNumMethod::setResult;

public:
    ScanningMethod(const OptProblemType &_problem = OptProblemType(),
                   const Parameters &parameters = Parameters()):
        opt::ICharacteristicOptMethod<Trial>(),
        GeneralNumMethod(_problem, parameters),
        cycleBoundary(0)
    {};

    void setParameters(const typename GeneralMethod::Parameters &parameters) override {
        GeneralNumMethod::setParameters(parameters);
    };
    void getParameters(typename GeneralMethod::Parameters &parameters) const override {
        GeneralNumMethod::getParameters(parameters);
    };

    void solve(typename GeneralMethod::Result &result) override;
    bool solveTest(typename GeneralMethod::Result &result) override;
};

template <typename OptProblemType>
void ScanningMethod<OptProblemType>::Report::printOptProblem(
    std::ostream &stream,
    const OptProblemType &optProblem) const
{
    stream << "[a; b] = [" << optProblem.getSearchArea().lowerBound << "; " <<
                              optProblem.getSearchArea().upBound << "]"<< "\n";
    std::vector<double> optimalPoints;
    optProblem.getOptimalPoints(optimalPoints);
    if (!optimalPoints.empty()) {
        stream << "X* = (" << optimalPoints[0];
        int numberOptimalPoints = optimalPoints.size();
        for (int i = 1; i < numberOptimalPoints; i++) {
            stream << "; " << optimalPoints[i];
        }
        stream << ")\n";
        stream << "f(X*) = " << optProblem.getOptimalValue() << "\n";
    }
}

template <typename OptProblemType>
void ScanningMethod<OptProblemType>::Report::printErrorEstimate(
    std::ostream &stream,
    const OptProblemType &optProblem,
    const typename GeneralMethod::Result &result) const
{
    double point = result.point;
    std::vector<double> optimalPoints;
    optProblem.getOptimalPoints(optimalPoints);
    if (!optimalPoints.empty()) {
        auto iter = std::min_element(optimalPoints.begin(), optimalPoints.end(),
            [&point] (const double &point1, const double &point2) {
                return std::abs(point1 - point) < std::abs(point2 - point);
            });
        stream << "|X* - X| = " << std::abs(*iter - point) << "\n";
        stream << "|f(X*) - f(X)| = " << std::abs(optProblem.getOptimalValue() - result.value) << "\n";
    }
}


template <typename OptProblemType>
void ScanningMethod<OptProblemType>::insertInSorted(const Trial &trial) {
    auto iter = this->trialPoints.begin();

    std::advance(iter, t);
    this->trialPoints.insert(iter, trial);
}

template <typename OptProblemType>
void ScanningMethod<OptProblemType>::calcCharacteristic() {
    t += 2;
    if (this->numberTrials - 2 == cycleBoundary) {
        t = 1;
        cycleBoundary <<= 1;
        cycleBoundary++; // TODO: check priority of operations
    }

    // Linear complexity
    // double R = -std::numeric_limits<double>::infinity(), Rtmp;
    // size_t trialPointsSize = this->trialPoints.size();
    // double xPrev = this->trialPoints[0].x, xNext;
    // for (size_t i = 1; i < trialPointsSize; i++) {
    //     xNext = this->trialPoints[i].x;
    //     Rtmp = xNext - xPrev;
    //     if (Rtmp > R) {
    //         R = Rtmp;
    //         t = (int)i;
    //     }
    //     xPrev = xNext;
    // }
}

template <typename OptProblemType>
Trial ScanningMethod<OptProblemType>::newTrial(const typename OptProblemType::Point &x) {
    this->numberFevals++;
    return Trial(x, this->problem.computeObjFunction(x));
}

template <typename OptProblemType>
typename OptProblemType::Point ScanningMethod<OptProblemType>::selectNewPoint() {
    calcCharacteristic();
    return 0.5 * (this->trialPoints[t].x + this->trialPoints[(size_t)t - 1].x);
}

template <typename OptProblemType>
double ScanningMethod<OptProblemType>::estimateSolution(typename OptProblemType::Point &x) const {
    auto iter = min_element(this->trialPoints.begin(), this->trialPoints.end(),
        [] (const Trial &trial1, const Trial &trial2) {
            return trial1.z < trial2.z;
        });
    x = iter->x;

    return iter->z;
}

template <typename OptProblemType>
bool ScanningMethod<OptProblemType>::stopConditions() {
    if (std::abs(this->trialPoints[t].x - this->trialPoints[(size_t)t - 1].x) <= this->accuracy) {
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
bool ScanningMethod<OptProblemType>::stopConditionsTest() {
    std::vector<double> optimalPoints;
    this->problem.getOptimalPoints(optimalPoints);
    int numberOptimalPoints = (int)optimalPoints.size();
    for (int i = 0; i < numberOptimalPoints; i++) {
        if (std::abs(this->trialPoints[t].x - optimalPoints[i]) <= this->error) {
            this->stoppingCondition = StoppingConditions::ERROR;
            return true;
        }
    }
    return stopConditions();
}

template <typename OptProblemType>
void ScanningMethod<OptProblemType>::solve(typename GeneralMethod::Result &result) {
    this->trialPoints.clear();
    this->numberFevals = 0;

    this->trialPoints.push_back(newTrial(this->problem.getSearchArea().lowerBound));
    this->trialPoints.push_back(newTrial(this->problem.getSearchArea().upBound));
    t = 1;

    this->numberTrials = 2;
    cycleBoundary = 0;

    Trial trial;
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
bool ScanningMethod<OptProblemType>::solveTest(typename GeneralMethod::Result &result) {
    this->trialPoints.clear();
    this->numberFevals = 0;

    this->trialPoints.push_back(newTrial(this->problem.getSearchArea().lowerBound));
    this->trialPoints.push_back(newTrial(this->problem.getSearchArea().upBound));
    t = 1;

    this->numberTrials = 2;
    cycleBoundary = 0;

    Trial trial;
    double xNew;
    while (stopConditionsTest()) {
        xNew = selectNewPoint();

        trial = newTrial(xNew);
        this->numberTrials++;

        insertInSorted(trial);
    }

    setResult(result);

    return static_cast<Result&>(result).stoppingCondition == StoppingConditions::ERROR;
}


#endif // SCANNING_METHOD_H_
