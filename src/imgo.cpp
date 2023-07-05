/* #include <imgo.h>

#include <algorithm>
#include <limits>
#if defined( _MSC_VER )
    #define _USE_MATH_DEFINES
    #include <math.h>
#else
    #include <cmath>
#endif

const double epsilon = 1e-14;

inline int insertInSorted(vector<TrialConstrained> &trials, TrialConstrained trial) {
    vector<TrialConstrained>::iterator iter = trials.begin();
    vector<TrialConstrained>::iterator iterEnd = trials.end();
    int pos = 0;
    while(true) {
        if (iter == iterEnd || iter->x > trial.x) break;
        iter++; pos++;
    }
    trials.insert(iter, trial);
    return pos;
}

inline double searchMin(vector<TrialConstrained> &trials, int numberConstraints) {
    double z = numeric_limits<double>::infinity(), x = 0.0;
    size_t sizeTrials = trials.size();
    for (int i = 0; i < sizeTrials; i++) {
        if (trials[i].nu == numberConstraints + 1 && trials[i].z < z) {
            z = trials[i].z;
            x = trials[i].x;
        }
    }
    return x;
}

TrialConstrained ImgoMethod::newTrial(double x) {
    TrialConstrained trial(x);
    for (int j = 1; j <= numberConstraints + 1; j++) {
        numberFevals++;
        if ((f(x, j) > 0) || (j == numberConstraints + 1)) {
            trial.z = f(x, j);
            trial.nu = j;
            break;
        }
    }
    return trial;
}

double ImgoMethod::newPoint(int t) {
    if (trialPoints[t].nu != trialPoints[(size_t)t - 1].nu) {
        return (trialPoints[t].x + trialPoints[(size_t)t - 1].x) / 2.0;
    } else {
        return (trialPoints[t].x + trialPoints[(size_t)t - 1].x) / 2.0 - 
               (trialPoints[t].z - trialPoints[(size_t)t - 1].z) / (2.0 * r * mu[(size_t)trialPoints[t].nu - 1]);
    }
}

double ImgoMethod::selectNewPoint(int &t) {
    int nuLastTrial;
    size_t sizeI;
    double muTmp;

    // Step 3
    // with optimization(const)
    nuLastTrial = lastTrial.nu - 1;
    sizeI = I[nuLastTrial].size();
    for (int nu = 0; nu < numberConstraints + 1; nu++) {
        if (!calcI[nu]) mu[nu] = 0.0;
    }
    if (I[nuLastTrial].size() >= 3) {
        if (lastTrialPos == 0) {
            mu[nuLastTrial] = max({ mu[nuLastTrial], abs(I[nuLastTrial][1].z - I[nuLastTrial][0].z) /
                                                     pow(I[nuLastTrial][1].x - I[nuLastTrial][0].x, 1.0 / n) });  
        } else if (lastTrialPos == I[nuLastTrial].size() - 1) {
            mu[nuLastTrial] = max({ mu[nuLastTrial],
                                    abs(I[nuLastTrial][sizeI - 1].z - I[nuLastTrial][sizeI - 2].z) / 
                                    pow(I[nuLastTrial][sizeI - 1].x - I[nuLastTrial][sizeI - 2].x, 1.0 / n) });
        } else {
            mu[nuLastTrial] = max({ mu[nuLastTrial],
                                    abs(I[nuLastTrial][lastTrialPos].z - I[nuLastTrial][(size_t)lastTrialPos - 1].z) / 
                                    pow(I[nuLastTrial][lastTrialPos].x - I[nuLastTrial][(size_t)lastTrialPos - 1].x, 1.0 / n),
                                    abs(I[nuLastTrial][(size_t)lastTrialPos + 1].z - I[nuLastTrial][lastTrialPos].z) / 
                                    pow(I[nuLastTrial][(size_t)lastTrialPos + 1].x - I[nuLastTrial][lastTrialPos].x, 1.0 / n) });
        }
    } else if (I[nuLastTrial].size() == 2) {
        mu[nuLastTrial] = max({ mu[nuLastTrial], abs(I[nuLastTrial][1].z - I[nuLastTrial][0].z) /
                                                 pow(I[nuLastTrial][1].x - I[nuLastTrial][0].x, 1.0 / n) });
    }
    if (abs(mu[nuLastTrial]) > epsilon) calcI[nuLastTrial] = true;
    for (int nu = 0; nu < numberConstraints + 1; nu++) {
        if (abs(mu[nu]) <= epsilon) mu[nu] = 1.0;
    }

    // with optimization(linear)
    // nuLastTrial = lastTrial.nu - 1;
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

    // Step 4
    for (int nu = 0; nu < numberConstraints + 1; nu++) {
        if (I[nu].size() != 0) {
            zStar[nu] = I[nu][0].z;
            sizeI = I[nu].size();
            for (int i = 1; i < sizeI; i++) {
                if (I[nu][i].z < zStar[nu]) {
                    zStar[nu] = I[nu][i].z;
                }
            }
            if (zStar[nu] <= 0.0 && nu != numberConstraints) {
                zStar[nu] = -mu[nu] * d;
            }
        }
    }

    // Steps 5, 6
    double R = -numeric_limits<double>::infinity(), Rtmp;
    double muV, zStarV, dx;

    size_t sizeTrialPoints = trialPoints.size();
    for (size_t i = 1; i < sizeTrialPoints; i++) {
        dx = trialPoints[i].x - trialPoints[i - 1].x;
        if (trialPoints[i].nu == trialPoints[i - 1].nu) {
            muV = mu[(size_t)trialPoints[i].nu - 1];
            zStarV = zStar[(size_t)trialPoints[i].nu - 1];
            Rtmp = dx + pow(trialPoints[i].z - trialPoints[i - 1].z, 2) / (r * r * muV * muV * dx) -
                   2.0 * (trialPoints[i].z + trialPoints[i - 1].z - 2.0 * zStarV) / (r * muV);
        } else if (trialPoints[i - 1].nu < trialPoints[i].nu) {
            muV = mu[(size_t)trialPoints[i].nu - 1];
            zStarV = zStar[(size_t)trialPoints[i].nu - 1];
            Rtmp = 2.0 * dx - 4.0 * (trialPoints[i].z - zStarV) / (r * muV);
        } else  {
            muV = mu[(size_t)trialPoints[i - 1].nu - 1];
            zStarV = zStar[(size_t)trialPoints[i - 1].nu - 1];
            Rtmp = 2.0 * dx - 4.0 * (trialPoints[i - 1].z - zStarV) / (r * muV);
        }
        if (Rtmp > R) {
            R = Rtmp;
            t = (int)i;
        }
    }

    // Step 7
    return newPoint(t);
}

void ImgoMethod::setNumberConstraints(int _numberConstraints) {
    OptimizationMethodConstrained::setNumberConstraints(_numberConstraints);
    I.resize((size_t)numberConstraints + 1);
    calcI.resize((size_t)numberConstraints + 1);
    mu.resize((size_t)numberConstraints + 1);
    zStar.resize((size_t)numberConstraints + 1);
}

void ImgoMethod::solve(int &numberTrials, int &numberFevals, double &x) {
    for (int i = 0; i < I.size(); i++) {
        I[i].clear();
        calcI[i] = false;
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
    while(true) {
        numberTrials++;

        // Steps 3, 4, 5, 6, 7
        xNew = selectNewPoint(t);
        lastTrial = newTrial(xNew);

        // Step 1
        insertInSorted(trialPoints, lastTrial);

        // Step 2
        lastTrialPos = insertInSorted(I[(size_t)lastTrial.nu - 1], lastTrial);

        // Stop conditions
        if (trialPoints[t].x - trialPoints[(size_t)t - 1].x <= eps) break;
        if (this->numberFevals >= maxFevals || numberTrials >= maxTrials) break;
    }
    numberFevals = this->numberFevals;
    x = searchMin(trialPoints, numberConstraints);
}

void ImgoMethod::solve(int &numberTrials, int &numberFevals, vector<double> &X) {
    solve(numberTrials, numberFevals, X[0]);
}

bool ImgoMethod::solveTest(double xOpt, int &numberTrials, int &numberFevals) {
    for (int nu = 0; nu < numberConstraints + 1; nu++) {
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
    }
}
 */