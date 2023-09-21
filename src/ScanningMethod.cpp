#include <opt_methods/ScanningMethod.h>

#include <algorithm>
#include <limits>
#include <iterator>

#include <my_math.h>

using std::numeric_limits;
using std::advance;
using std::distance;
using std::min_element;
using std::max;

const double epsilon = 1e-14;

void ScanningMethod::Report::printOptProblem(
    std::ostream &stream,
    const OneDimensionalProblem &optProblem) const
{
    stream << "[a; b] = [" << optProblem.getSearchArea().lowerBound << "; " <<
                              optProblem.getSearchArea().upBound << "]"<< "\n";
    vector<double> optimalPoints;
    optProblem.getOptimalPoints(optimalPoints);
    stream << "X* = (" << optimalPoints[0];
    int numberOptimalPoints = optimalPoints.size();
    for (int i = 1; i < numberOptimalPoints; i++) {
        stream << "; " << optimalPoints[i];
    }
    stream << ")\n";
    stream << "f(X*) = " << optProblem.getOptimalValue() << "\n";
}

void ScanningMethod::Report::printErrorEstimate(
    std::ostream &stream,
    const OneDimensionalProblem &optProblem,
    const GeneralMethod::Result &result) const
{
    double point = result.point;
    vector<double> optimalPoints;
    optProblem.getOptimalPoints(optimalPoints);
    auto iter = std::min_element(optimalPoints.begin(), optimalPoints.end(),
        [&point] (const double &point1, const double &point2) {
            return std::abs(point1 - point) < std::abs(point2 - point);
        });
    stream << "|X* - X| = " << std::abs(*iter - point) << "\n";
    stream << "|f(X*) - f(X)| = " << std::abs(optProblem.getOptimalValue() - result.value) << "\n";
}


void ScanningMethod::insertInSorted(const Trial &trial) {
    vector<Trial>::iterator iter = trialPoints.begin();

    advance(iter, t);
    trialPoints.insert(iter, trial);
}

void ScanningMethod::calcCharacteristic() {
    t += 2;
    if (numberTrials - 2 == cycleBoundary) {
        t = 1;
        cycleBoundary <<= 1;
        cycleBoundary++; // TODO: check priority of operations
    }

    // Linear complexity
    // double R = -numeric_limits<double>::infinity(), Rtmp;
    // size_t trialPointsSize = trialPoints.size();
    // double xPrev = trialPoints[0].x, xNext;
    // for (size_t i = 1; i < trialPointsSize; i++) {
    //     xNext = trialPoints[i].x;
    //     Rtmp = xNext - xPrev;
    //     if (Rtmp > R) {
    //         R = Rtmp;
    //         t = (int)i;
    //     }
    //     xPrev = xNext;
    // }
}

Trial ScanningMethod::newTrial(const OneDimensionalProblem::Point &x) {
    numberFevals++;
    return Trial(x, problem.computeObjFunction(x));
}

double ScanningMethod::selectNewPoint() {
    calcCharacteristic();

    return 0.5 * (trialPoints[t].x + trialPoints[(size_t)t - 1].x);
}

double ScanningMethod::estimateSolution(OneDimensionalProblem::Point &x) const {
    auto iter = min_element(trialPoints.begin(), trialPoints.end(),
        [] (const Trial &trial1, const Trial &trial2) {
            return trial1.z < trial2.z;
        });
    x = iter->x;

    return iter->z;
}

bool ScanningMethod::stopConditions() {
    if (std::abs(trialPoints[t].x - trialPoints[(size_t)t - 1].x) <= accuracy) {
        stoppingCondition = StoppingConditions::ACCURACY;
        return true;
    }
    if (numberTrials >= maxTrials) {
        stoppingCondition = StoppingConditions::MAXTRIALS;
        return true;
    }
    if (numberFevals >= maxFevals) {
        stoppingCondition = StoppingConditions::MAXFEVALS;
        return true;
    }
    return false;
}

bool ScanningMethod::stopConditionsTest() {
    vector<double> optimalPoints;
    problem.getOptimalPoints(optimalPoints);
    int numberOptimalPoints = (int)optimalPoints.size();
    for (int i = 0; i < numberOptimalPoints; i++) {
        if (std::abs(trialPoints[t].x - optimalPoints[i]) <= error) {
            stoppingCondition = StoppingConditions::ERROR;
            return true;
        }
    }
    return stopConditions();
}

void ScanningMethod::solve(GeneralMethod::Result &result) {
    trialPoints.clear();
    numberFevals = 0;

    trialPoints.push_back(newTrial(problem.getSearchArea().lowerBound));
    trialPoints.push_back(newTrial(problem.getSearchArea().upBound));
    t = 1;

    numberTrials = 2;
    cycleBoundary = 0;

    Trial trial;
    double xNew;
    while(!stopConditions()) {
        xNew = selectNewPoint();

        trial = newTrial(xNew);
        numberTrials++;

        insertInSorted(trial);
    }

    setResult(result);
}

bool ScanningMethod::solveTest(GeneralMethod::Result &result) {
    trialPoints.clear();
    numberFevals = 0;

    trialPoints.push_back(newTrial(problem.getSearchArea().lowerBound));
    trialPoints.push_back(newTrial(problem.getSearchArea().upBound));
    t = 1;

    numberTrials = 2;
    cycleBoundary = 0;

    Trial trial;
    double xNew;
    while (stopConditionsTest()) {
        xNew = selectNewPoint();

        trial = newTrial(xNew);
        numberTrials++;

        insertInSorted(trial);
    }

    setResult(result);

    return static_cast<Result&>(result).stoppingCondition == StoppingConditions::ERROR;
}
