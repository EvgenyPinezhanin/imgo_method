#include <opt_methods/GsaMethod.h>

#include <algorithm>
#include <limits>
#include <iterator>

#include <result_methods/StronginResultMethod.h>
#include <MyMath.h>

using std::numeric_limits;
using std::advance;
using std::distance;
using std::max;
using std::abs;

const double epsilon = 1e-14;

int GsaMethod::insertInSorted(Trial trial) {
    vector<Trial>::iterator iter = trialPoints.begin();
    int dist = trialPoints.size();

    advance(iter, dist / 2);
    while (true) {
        dist /= 2;
        if (trial.x < iter->x) advance(iter, -dist / 2);
        else {
            if (trial.x < (iter + 1)->x) break;
            else advance(iter, dist / 2);
        }
    }
    iter = trialPoints.insert(iter, trial);

    return distance(trialPoints.begin(), iter);
}

double GsaMethod::searchMin() const {
    double z = numeric_limits<double>::infinity(), x = 0.0;
    size_t trialPointsSize = trialPoints.size();
    for (int i = 0; i < trialPointsSize; i++) {
        if (trialPoints[i].z < z) {
            z = trialPoints[i].z;
            x = trialPoints[i].x;
        }
    }
    return x;
}

void GsaMethod::calcCharacteristic() {
    double R = -numeric_limits<double>::infinity(), Rtmp, dx;
    
    int sizeTrialPoints = trialPoints.size();
    for (size_t i = 1; i < sizeTrialPoints; i++) {
        dx = trialPoints[i].x - trialPoints[i - 1].x;
        Rtmp = constantEstimation * dx + pow(trialPoints[i].z - trialPoints[i - 1].z, 2) / (constantEstimation * dx) -
               2 * (trialPoints[i].z + trialPoints[i - 1].z);
        if (Rtmp > R) {
            R = Rtmp;
            t = (int)i;
        }
    }
}

Trial GsaMethod::newTrial(const double &x) {
    numberFevals++;
    return Trial(x, task.computeObjFunction(x));
}

double GsaMethod::newPoint() {
    return (trialPoints[t].x + trialPoints[(size_t)t - 1].x) / 2 - (trialPoints[t].z - trialPoints[(size_t)t - 1].z) / (2 * constantEstimation);
}

double GsaMethod::selectNewPoint() {
    static double M = -1.0;
    size_t sizeTrialPoints;

    // Step 2
    // with optimization(const)
    if (lastTrials[0].x == task.getSearchArea().getUpBound()) {
        M = abs((trialPoints[1].z - trialPoints[0].z) / (trialPoints[1].x - trialPoints[0].x));
    } else {
        M = max({ M, abs((lastTrials[0].z - trialPoints[(size_t)lastTrialsPos[0] - 1].z) / 
                         (lastTrials[0].x - trialPoints[(size_t)lastTrialsPos[0] - 1].x)), 
                     abs((trialPoints[(size_t)lastTrialsPos[0] + 1].z - lastTrials[0].z) / 
                         (trialPoints[(size_t)lastTrialsPos[0] + 1].x - lastTrials[0].x)) });
    }

    // without optimization
    // double Mtmp;
    // sizeTrialPoints = trialPoints.size();
    // M = 0.0;
    // for (int i = 1; i < sizeTrialPoints; i++) {
    //     Mtmp = abs(trialPoints[i].z - trialPoints[(size_t)i - 1].z) / (trialPoints[i].x - trialPoints[(size_t)i - 1].x);
    //     if (MTmp > M) M = Mtmp;
    // }

    // Step 3
    constantEstimation = (abs(M) <= epsilon) ? 1.0 : r * M;

    // Steps 4, 5
    calcCharacteristic();

    // Step 6
    return newPoint();
}

inline double GsaMethod::estimateSolution() const {
    return searchMin();
}

bool GsaMethod::stopConditions() {
    if (trialPoints[t].x - trialPoints[(size_t)t - 1].x <= accuracy) {
        stopCriteria = StopCriteria::accuracy;
        return true;
    }
    if (numberFevals >= maxFevals) {
        stopCriteria = StopCriteria::maxFevals;
        return true;
    }
    if (numberTrials >= maxTrials) {
        stopCriteria = StopCriteria::maxTrials;
        return true;
    }
    return false;
}

void GsaMethod::solve(StronginResultMethod<double> &result) {
    trialPoints.clear();
    numberFevals = 0;

    trialPoints.push_back(newTrial(task.getSearchArea().getLowerBound()));
    lastTrials[0] = newTrial(task.getSearchArea().getUpBound());

    numberTrials = 2;

    double xNew;
    while(true) {
        // Step 1
        lastTrialsPos[0] = insertInSorted(lastTrials[0]);

        // Steps 2, 3, 4, 5, 6
        xNew = selectNewPoint();

        lastTrials[0] = newTrial(xNew);
        numberTrials++;

        if (stopConditions()) break;
    }
    setDataInResultMethod(result);
}

bool GsaMethod::solveTest(StronginResultMethod<double> &result) {
/*     for (int nu = 0; nu < numberConstraints + 1; nu++) {
        I[nu].clear();
        calcI[nu] = false;
    }
    trialPoints.clear();
    this->numberFevals = 0;

    lastTrial = newTrial(A[0]);
    trialPoints.push_back(lastTrial);
    insertInSorted(I[(size_t)lastTrial.nu - 1], lastTrial);
    lastTrial = newTrial(B[0]);
    trialPoints.push_back(lastTrial);
    lastTrialPos = insertInSorted(I[(size_t)lastTrial.nu - 1], lastTrial);
    numberTrials = 2;

    double xNew;
    int t;
    while (true) {
        numberTrials++;

        // Steps 3, 4, 5, 6, 7
        xNew = selectNewPoint(t);
        lastTrial = newTrial(xNew);

        // Step 1
        insertInSorted(trialPoints, lastTrial);

        // Step 2
        lastTrialPos = insertInSorted(I[(size_t)lastTrial.nu - 1], lastTrial);

        // Stop conditions
        if (abs(xNew - xOpt) <= eps) {
            numberFevals = this->numberFevals;
            return true;
        }
        if (this->numberFevals >= maxFevals || numberTrials >= maxTrials) {
            numberFevals = this->numberFevals;
            return false;
        }
    } */
    return false;
}
