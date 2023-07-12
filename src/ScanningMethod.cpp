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
using std::abs;

const double epsilon = 1e-14;

void ScanningMethod::insertInSorted(const Trial &trial) {
    vector<Trial>::iterator iter = trialPoints.begin();

    advance(iter, t);
    trialPoints.insert(iter, trial);
}

void ScanningMethod::calcCharacteristic() {
    double R = -numeric_limits<double>::infinity(), Rtmp;
    
    int trialPointsSize = (int)trialPoints.size();
    for (int i = 1; i < trialPointsSize; i++) {
        Rtmp = trialPoints[i].x - trialPoints[i - 1].x;
        if (Rtmp > R) {
            R = Rtmp;
            t = i;
        }
    }
}

Trial ScanningMethod::newTrial(const double &x) {
    numberFevals++;
    return Trial(x, task.computeObjFunction(x));
}

double ScanningMethod::selectNewPoint() {
    // Steps 2, 3
    calcCharacteristic();

    // Step 4
    return 0.5 * (trialPoints[t].x + trialPoints[(size_t)t - 1].x);
}

double ScanningMethod::estimateSolution() const {
    auto iter = min_element(trialPoints.begin(), trialPoints.end(),
    [] (const Trial &trial1, const Trial &trial2) {
        return trial1.z < trial2.z;
    });

    return iter->x;
}

bool ScanningMethod::stopConditions() {
    if (abs(trialPoints[t].x - trialPoints[(size_t)t - 1].x) <= accuracy) {
        result.stoppingCondition = StoppingCondition::accuracy;
        return true;
    }
    if (numberFevals >= maxFevals) {
        result.stoppingCondition = StoppingCondition::maxFevals;
        return true;
    }
    if (numberTrials >= maxTrials) {
        result.stoppingCondition = StoppingCondition::maxTrials;
        return true;
    }
    return false;
}

bool ScanningMethod::stopConditionsTest() {
    if (abs(trialPoints[t].x - task.getOptimalPoint()) <= error) {
        result.stoppingCondition = StoppingCondition::error;
        return true;
    }
    return stopConditions();
}

void ScanningMethod::solve(ResultMethod &_result) {
    trialPoints.clear();
    numberFevals = 0;

    // Step 1
    trialPoints.push_back(newTrial(task.getSearchArea().getLowerBound()));
    trialPoints.push_back(newTrial(task.getSearchArea().getUpBound()));
    t = 1;

    numberTrials = 2;

    Trial trial;
    double xNew;
    while(!stopConditions()) {
        // Steps 2, 3, 4
        xNew = selectNewPoint();

        trial = newTrial(xNew);
        numberTrials++;

        // Step 1
        insertInSorted(trial);
    }

    result.numberTrials = numberTrials;
    result.numberFevals = numberFevals;
    result.solution = estimateSolution();
    _result = result;
}

bool ScanningMethod::solveTest(ResultMethod &_result) {
    trialPoints.clear();
    numberFevals = 0;

    result = _result;

    trialPoints.push_back(newTrial(task.getSearchArea().getLowerBound()));
    Trial trial = newTrial(task.getSearchArea().getUpBound());

    numberTrials = 2;

    double xNew;
    while (true) {
        // Step 1
        insertInSorted(trial);

        // Steps 2, 3, 4
        xNew = selectNewPoint();

        trial = newTrial(xNew);
        numberTrials++;

        if (stopConditionsTest()) break;
    }
    
    result.numberTrials = numberTrials;
    result.numberFevals = numberFevals;
    result.solution = estimateSolution();

    return result.stoppingCondition == StoppingCondition::error;
}
