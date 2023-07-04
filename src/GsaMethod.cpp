#include <GsaMethod.h>

#include <algorithm>
#include <limits>
#include <iterator>

#include <StronginResultMethod.h>
#include <MyMath.h>

using std::numeric_limits;
using std::advance;
using std::distance;

const double epsilon = 1e-14;

double GsaMethod::compute(vector<double> X) const {
    return objFunction(X[0]);
}

int GsaMethod::insertInSorted(vector<Trial> &trials, Trial trial) {
    vector<Trial>::iterator iter = trials.begin();
    int dist = trials.size();

    advance(iter, dist / 2);
    while (true) {
        dist /= 2;
        if (trial.x < iter->x) advance(iter, -dist / 2);
        else {
            if (trial.x < (iter + 1)->x) break;
            else advance(iter, dist / 2);
        }
    }
    iter = trials.insert(iter, trial);

    return distance(trials.begin(), iter);
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

Trial GsaMethod::newTrial(double x) {
    numberFevals++;
    return Trial(x, f(x));
}

double GsaMethod::newPoint(int t) {
    return (trialPoints[t].x + trialPoints[(size_t)t - 1].x) / 2 - (trialPoints[t].z - trialPoints[(size_t)t - 1].z) / (2 * m);
}

double GsaMethod::selectNewPoint(int &t) {
    static double M = -1.0;
    size_t sizeTrialPoints;

    // Step 2
    // with optimization(const)
    if (lastTrial.x == B[0]) {
        M = abs((trialPoints[1].z - trialPoints[0].z) / (trialPoints[1].x - trialPoints[0].x));
    } else {
        M = max({ M, abs((lastTrial.z - trialPoints[(size_t)lastTrialPos - 1].z) / 
                         (lastTrial.x - trialPoints[(size_t)lastTrialPos - 1].x)), 
                     abs((trialPoints[(size_t)lastTrialPos + 1].z - lastTrial.z) / 
                         (trialPoints[(size_t)lastTrialPos + 1].x - lastTrial.x)) });
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
    m = (abs(M) <= epsilon) ? 1.0 : r * M;

    // Steps 4, 5
    double R = -numeric_limits<double>::infinity(), Rtmp, dx;
    
    sizeTrialPoints = trialPoints.size();
    for (size_t i = 1; i < sizeTrialPoints; i++) {
        dx = trialPoints[i].x - trialPoints[i - 1].x;
        Rtmp = m * dx + pow(trialPoints[i].z - trialPoints[i - 1].z, 2) / (m * dx) - 2 * (trialPoints[i].z + trialPoints[i - 1].z);
        if (Rtmp > R) {
            R = Rtmp;
            t = (int)i;
        }
    }

    // Step 6
    return newPoint(t);
}

void GsaMethod::solve(int &numberIters, int &numberFevals, double &x) {
    trialPoints.clear();
    this->numberFevals = 0;

    trialPoints.push_back(newTrial(A[0]));

    lastTrial = newTrial(B[0]);
    lastTrialPos = 1;

    trialPoints.push_back(lastTrial);
    numberIters = 2;

    double xNew;
    int t;
    while(true) {
        numberIters++;

        // Steps 2, 3, 4, 5, 6
        xNew = selectNewPoint(t);

        // Trial
        lastTrial = newTrial(xNew);

        // Stop conditions
        if (trialPoints[t].x - trialPoints[(size_t)t - 1].x <= eps) break;
        if (this->numberFevals >= maxFevals || numberIters >= maxTrials) break;

        // Step 1
        lastTrialPos = insertInSorted(trialPoints, lastTrial);
    }
    numberFevals = this->numberFevals;
    x = searchMin(trialPoints);
}

void GsaMethod::solve(int &numberTrials, int &numberFevals, vector<double> &X) {
    solve(numberTrials, numberFevals, X[0]);
}
