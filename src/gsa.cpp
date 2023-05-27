#include <gsa.h>

#include <algorithm>
#include <limits>
#if defined( _MSC_VER )
    #define _USE_MATH_DEFINES
    #include <math.h>
#else
    #include <cmath>
#endif

const double epsilon = 1e-14;

inline int insertInSorted(vector<Trial> &trials, Trial trial) {
    vector<Trial>::iterator iter = trials.begin();
    vector<Trial>::iterator iterEnd = trials.end();
    int pos = 0;
    while(true) {
        if (iter == iterEnd || iter->x > trial.x) break;
        iter++; pos++;
    }
    trials.insert(iter, trial);
    return pos;
}

inline double searchMin(vector<Trial> &trials) {
    double z = numeric_limits<double>::infinity(), x = 0.0;
    size_t sizeTrials = trials.size();
    for (int i = 0; i < sizeTrials; i++) {
        if (trials[i].z < z) {
            z = trials[i].z;
            x = trials[i].x;
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
