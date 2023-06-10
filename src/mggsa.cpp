#include <mggsa.h>

#include <iostream>
#include <fstream>
#include <ctime>
#include <limits>
#include <algorithm>
#if defined( _MSC_VER )
    #define _USE_MATH_DEFINES
    #include <math.h>
#else
    #include <cmath>
#endif

#include <map.h>
#include <output_results.h>
#include <omp.h>

// #define DEBUG
// #define EPS
// #define TIME_TEST

const double epsilon = 1e-14;

#if defined(TIME_TEST)
    ofstream ofstrTimeTest;
    vector<int> numberTimeStamp;
    vector<double> muTime, zStarTime, RTime, trialTime, addTrialTime, addITime, PTime;
    double startTime, endTime;
#endif

const double peanoA = 0.0, peanoB = 1.0, peanoRandom = 0.5;

inline double euclideanDistance(vector<double> val1, vector<double> val2) {
    double res = 0.0;
    size_t size = val1.size();
    for (int i = 0; i < size; i++) {
        res += (val1[i] - val2[i]) * (val1[i] - val2[i]);
    }
    return sqrt(res);
}

template<typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

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

TrialConstrained MggsaMethod::newTrial(double x) {
    TrialConstrained trial(x);
    vector<double> X(n);
    y(x, X);
    for (int j = 1; j <= numberConstraints + 1; j++) {
        numberFevals++;
        if ((f(X, j) > 0) || (j == numberConstraints + 1)) {
            trial.z = f(X, j);
            trial.nu = j;
            break;
        }
    }

#if defined(DEBUG)
    cout << "trial: x = " << x << " X = " << X[0] << " Y = " << X[1] <<
            " z = " << f(X, numberConstraints + 1) << " nu = " << trial.nu << endl;
#endif

    return trial;
}

double MggsaMethod::newPoint(int t) {
#if defined(DEBUG)
    cout << "t = " << t << endl;
    cout << "trialPoints[t].x = " << trialPoints[t].x << " trialPoints[t-1].x = " << trialPoints[(size_t)t - 1].x << endl;
    cout << "mu[trialPoints[t].nu] = " << mu[(size_t)trialPoints[t].nu - 1] << endl;
#endif

    if (trialPoints[t].nu != trialPoints[(size_t)t - 1].nu) {
        return (trialPoints[t].x + trialPoints[(size_t)t - 1].x) / 2.0;
    } else {
        return (trialPoints[t].x + trialPoints[(size_t)t - 1].x) / 2.0 -
               sgn(trialPoints[t].z - trialPoints[(size_t)t - 1].z) / (2.0 * r) *
               pow(abs(trialPoints[t].z - trialPoints[(size_t)t - 1].z) / mu[(size_t)trialPoints[t].nu - 1], n);
    }
}

double MggsaMethod::calcH(double xNew) {
    double d = 1.0 / (pow(2.0, den * n) * (pow(2.0, n) - 1.0));
    double h = floor(xNew / d) * d;
    return h;
}

bool MggsaMethod::checkDensity(double h) {
    size_t sizeTrialPoints = trialPoints.size();
    for (int i = 0; i < sizeTrialPoints; i++) {
        if (abs(h - trialPoints[i].x) <= epsilon) return false; 
    }
    return true;
}

void MggsaMethod::y(double x, vector<double> &X) {
    int d = (key != 3) ? den : den + 1;
    mapd(x, d, X.data(), n, key);

#if defined(DEBUG)
    cout << "X = (";
    for (int i = 0; i < n - 1; i++) {
        cout << X[i] << ", ";
    }
    cout << X[n - 1] << ")" << endl;
#endif

    for (int i = 0; i < n; i++) {
        X[i] = X[i] * (B[i] - A[i]) + (A[i] + B[i]) / 2.0;
    }
}

void MggsaMethod::x(const vector<double> &P, vector<double> &X) {
    vector<double> pCorrect(n);
    int sizeX;
    X.resize((size_t)pow(2, n));
    for (int i = 0; i < n; i++) {
        pCorrect[i] = (P[i] - (A[i] + B[i]) / 2.0) / (B[i] - A[i]);
    }

    invmad(den + 1, X.data(), (int)pow(2, n), &sizeX, pCorrect.data(), n, incr);
    X.resize(sizeX);
}

double MggsaMethod::selectNewPoint(int &t) {
#if defined(TIME_TEST)
    startTime = omp_get_wtime();
#endif

    double muTmp;
    int nuLastTrials;
    size_t sizeI, sizeLastTrials;

    // Step 3
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

    // with optimization (linear)
    sizeLastTrials = lastTrials.size();
    nuLastTrials = lastTrials[0].nu - 1;
    sizeI = I[nuLastTrials].size();
    for (int nu = 0; nu < numberConstraints + 1; nu++) {
        if (!calcI[nu]) mu[nu] = 0.0;
    }
    for (int i = 0; i < sizeLastTrials; i++) {
        for (int j = 0; j < sizeI; j++) {
            if (I[nuLastTrials][j].x != lastTrials[i].x) {
                muTmp = abs(I[nuLastTrials][j].z - lastTrials[i].z) /
                        pow(abs(I[nuLastTrials][j].x - lastTrials[i].x), 1.0 / n);
                if (muTmp > mu[nuLastTrials]) {
                    mu[nuLastTrials] = muTmp;
                    if (abs(mu[nuLastTrials]) > epsilon) calcI[nuLastTrials] = true;
                }
            }
        }
    }
    for (int nu = 0; nu < numberConstraints + 1; nu++) {
        if (abs(mu[nu]) <= epsilon) mu[nu] = 1.0;
    }

    // without optimization
    // for (int nu = 0; nu < numberConstraints + 1; nu++) {
    //     mu[nu] = 0.0;
    // }
    // for (int nu = 0; nu < numberConstraints + 1; nu++) {
    //     sizeI = I[nu].size();
    //     for (int i = 1; i < sizeI; i++) { // при i = 0 - нет j
    //         for (int j = 0; j < i; j++) {
    //             muTmp = abs(I[nu][i].z - I[nu][j].z) / pow(I[nu][i].x - I[nu][j].x, 1.0 / n);
    //             if (muTmp > mu[nu]) {
    //                 mu[nu] = muTmp;
    //             }
    //         }
    //     }
    //     if (abs(mu[nu]) <= epsilon) {
    //         mu[nu] = 1.0;
    //     };
    // }

#if defined(TIME_TEST)
    endTime = omp_get_wtime();
    muTime.push_back(endTime - startTime);
#endif

#if defined(DEBUG)
    cout << "mu: ";
    for (int nu = 0; nu < numberConstraints + 1; nu++) {
        cout << mu[nu] << " ";
    }
    cout << endl;
#endif

#if defined(TIME_TEST)
    startTime = omp_get_wtime();
#endif

    // Step 4
    for (int nu = 0; nu < numberConstraints + 1; nu++) {
        if (I[nu].size() != 0) {
            if (nu + 1 != M) {
                zStar[nu] = -mu[nu] * d;
            } else {
                zStar[nu] = I[nu][0].z;
                sizeI = I[nu].size();
                for (int i = 1; i < sizeI; i++) {
                    if (I[nu][i].z < zStar[nu]) {
                        zStar[nu] = I[nu][i].z;
                    }
                }
            }
        }
    }

#if defined(TIME_TEST)
    endTime = omp_get_wtime();
    zStarTime.push_back(endTime - startTime);
#endif

#if defined(DEBUG)
    cout << "zStar: ";
    for (int nu = 0; nu < numberConstraints + 1; nu++) {
        cout << zStar[nu] << " ";
    }
    cout << endl;
#endif

#if defined(TIME_TEST)
    startTime = omp_get_wtime();
#endif

    // Steps 5, 6
    double R = -numeric_limits<double>::infinity(), Rtmp;
    double muV, zStarV, dx;

    size_t sizeTrialPoints = trialPoints.size();
    for (size_t i = 1; i < sizeTrialPoints; i++) {
        dx = pow(trialPoints[i].x - trialPoints[i - 1].x, 1.0 / n);
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

#if defined(DEBUG)
    cout << "R[" << trialPoints[i - 1].x << ", " << trialPoints[i].x << "] = " << Rtmp << endl;
#endif

    }

#if defined(TIME_TEST)
    endTime = omp_get_wtime();
    RTime.push_back(endTime - startTime);
#endif

    // Pre-Step 7
    return newPoint(t);
}

void MggsaMethod::setNumberConstraints(int _numberConstraints) {
    OptimizationMethodConstrained::setNumberConstraints(_numberConstraints);
    I.resize((size_t)numberConstraints + 1);
    calcI.resize((size_t)numberConstraints + 1);
    mu.resize((size_t)numberConstraints + 1);
    zStar.resize((size_t)numberConstraints + 1);
}

void MggsaMethod::solve(int &numberTrials, int &numberFevals, vector<double> &X, TypeSolve type) {
#if defined(TIME_TEST)
    ofstrTimeTest.open("output_data/mggsa_test_time.txt");
    if (!ofstrTimeTest.is_open()) cerr << "File opening error\n";
    if (type == TypeSolve::SOLVE) {
        numberTimeStamp.clear();
        muTime.clear();
        zStarTime.clear();
        RTime.clear();
        trialTime.clear();
        addTrialTime.clear();
        addITime.clear();
        PTime.clear();
    }
#endif

    if (type == TypeSolve::SOLVE) {
        for (int nu = 0; nu < numberConstraints + 1; nu++) {
            I[nu].clear();
            calcI[nu] = false;
        }
        points.clear();
        trialPoints.clear();
        lastTrials.resize(1);
        lastTrialsPos.resize(1);
        this->numberFevals = 0;
    }
    X.resize(n);
    coincideX = false;

    if (type == TypeSolve::SOLVE) {
        lastTrials[0] = TrialConstrained{peanoA, -1.0, 0};
        trialPoints.push_back(lastTrials[0]);
        lastTrials[0].x = peanoB;
        trialPoints.push_back(lastTrials[0]);
    }

#if defined(DEBUG)
    if (type == TypeSolve::SOLVE)
        cout << "Trial number: " << 1 << endl;
#endif

#if defined(TIME_TEST)
    startTime = omp_get_wtime();
#endif

    if (type == TypeSolve::SOLVE) {
        lastTrials[0] = newTrial(peanoRandom);
        numberTrials = 1;
    }

#if defined(TIME_TEST)
    endTime = omp_get_wtime();
#endif

#if defined(TIME_TEST)
    if (type == TypeSolve::SOLVE) {
        numberTimeStamp.push_back(0);
        muTime.push_back(0.0);
        zStarTime.push_back(0.0);
        RTime.push_back(0.0);
        PTime.push_back(0.0);
        trialTime.push_back(endTime - startTime);
    }
#endif

#if defined(TIME_TEST)
    startTime = omp_get_wtime();
#endif

    if (type == TypeSolve::SOLVE) {
        insertInSorted(trialPoints, lastTrials[0]);
    }

#if defined(TIME_TEST)
    endTime = omp_get_wtime();
    if (type == TypeSolve::SOLVE)
        addTrialTime.push_back(endTime - startTime);
#endif

#if defined(TIME_TEST)
    startTime = omp_get_wtime();
#endif

    if (type == TypeSolve::SOLVE) {
        lastTrialsPos[0] = insertInSorted(I[(size_t)lastTrials[0].nu - 1], lastTrials[0]);
    }

#if defined(TIME_TEST)
    endTime = omp_get_wtime();
    if (type == TypeSolve::SOLVE)
        addITime.push_back(endTime - startTime);
#endif

    if (type == TypeSolve::SOLVE) {
        M = lastTrials[0].nu;
    }

    double xNew, h, deltaT;
    int t;
    TrialConstrained trial;
    vector<double> P(n), point;
    while(true) {
    #if defined(TIME_TEST)
        numberTimeStamp.push_back(numberTrials);
    #endif

    #if defined(DEBUG)
        cout << "Trial number: " << numberTrials << endl;
    #endif

        // Steps 3, 4, 5, 6, Pre-7
        xNew = selectNewPoint(t);

    #if defined(DEBUG)
        cout << "x^{k+1} = " << xNew << endl;
    #endif

        if (key == 3) {
        #if defined(TIME_TEST)
            startTime = omp_get_wtime();
        #endif

            h = calcH(xNew);
            y(h, P);

        #if defined(DEBUG)
            cout << "h = " << h << endl;
            cout << "P = (";
            for (int i = 0; i < n - 1; i++) {
                cout << P[i] << ", ";
            }
            cout << P[n - 1] << ")" << endl;
        #endif

        #if defined(TIME_TEST)
            endTime = omp_get_wtime();
            PTime.push_back(endTime - startTime);
        #endif
        }

        deltaT = pow(trialPoints[t].x - trialPoints[(size_t)t - 1].x, 1.0 / n);

    #if defined(TIME_TEST)
        startTime = omp_get_wtime();
    #endif

        // Step 7
        if (key != 3) {
            lastTrials[0] = newTrial(xNew);
        } else {
            x(P, hNu);
            if (!checkDensity(hNu[0])) {
                coincideX = true;
                break;
            }
            lastTrials.clear();
            trial = newTrial(hNu[0]);
            for (int i = 0; i < hNu.size(); i++) {
                trial.x = hNu[i];
                lastTrials.push_back(trial);
            }

        #if defined(DEBUG)
            for (int i = 0; i < hNu.size(); i++) {
                cout << "hNu[" << i << "] = " << hNu[i] << " ";
            }
            cout << endl;
        #endif

        }

        numberTrials++;

    #if defined(TIME_TEST)
        endTime = omp_get_wtime();
        trialTime.push_back(endTime - startTime);
    #endif

    if (key != 3) {
        y(xNew, P);
    }
    point = P;
    point.push_back(lastTrials[0].z);
    points.push_back(point);

    #if defined(TIME_TEST)
        startTime = omp_get_wtime();
    #endif

        // Step 1
        for (int i = 0; i < lastTrials.size(); i++) {
            insertInSorted(trialPoints, lastTrials[i]);
        }

    #if defined(TIME_TEST)
        endTime = omp_get_wtime();
        addTrialTime.push_back(endTime - startTime);
    #endif

    #if defined(TIME_TEST)
        startTime = omp_get_wtime();
    #endif

        // Step 2
        if (key != 3) {
            lastTrialsPos[0] = insertInSorted(I[(size_t)lastTrials[0].nu - 1], lastTrials[0]);
        } else {
            lastTrialsPos.resize(lastTrials.size());
            for (int i = 0; i < lastTrials.size(); i++) {
                lastTrialsPos[i] = insertInSorted(I[(size_t)lastTrials[0].nu - 1], lastTrials[i]);
            }
        }

    #if defined(TIME_TEST)
        endTime = omp_get_wtime();
        addITime.push_back(endTime - startTime);
    #endif

        if (lastTrials[0].nu > M) {
            M = lastTrials[0].nu;
        }

    #if defined(EPS)
        cout << pow(abs(trialPoints[t].x - trialPoints[t - 1].x), 1.0 / n) << endl;
    #endif

        // Stop conditions
        if (deltaT <= eps) break;
        if (this->numberFevals >= maxFevals || numberTrials >= maxTrials) break;
    }
    numberFevals = this->numberFevals;
    y(searchMin(trialPoints, numberConstraints), X);

#if defined(TIME_TEST)
    int k = (key == 3);
    for (int i = 0; i < numberTrials; i++) {
        ofstrTimeTest << numberTimeStamp[i] << " " << muTime[i] << " " << zStarTime[i] << " "
                      << RTime[i] << " " << trialTime[i] << " "  << addTrialTime[i] << " " << addITime[i];
        if (k) ofstrTimeTest << " " << PTime[i];
        ofstrTimeTest << endl;
    }
    ofstrTimeTest.close();

    drawGraphGnuplot("scripts/mggsa_test_time.gp", k);
#endif
}

void MggsaMethod::solve(int &numberTrials, int &numberFevals, vector<double> &X) {
    solve(numberTrials, numberFevals, X, TypeSolve::SOLVE);
}

bool MggsaMethod::solveTest(vector<double> XOpt, int &numberTrials, int &numberFevals, TypeSolve type) {
    if (type == TypeSolve::SOLVE) {
        for (int nu = 0; nu < numberConstraints + 1; nu++) {
            I[nu].clear();
            calcI[nu] = false;
        }
        trialPoints.clear();
        lastTrials.resize(1);
        lastTrialsPos.resize(1);
        this->numberFevals = 0;

        lastTrials[0] = TrialConstrained{peanoA, -1.0, 0};
        trialPoints.push_back(lastTrials[0]);
        lastTrials[0].x = peanoB;
        trialPoints.push_back(lastTrials[0]);

        lastTrials[0] = newTrial(peanoRandom);
        numberTrials = 1;

        insertInSorted(trialPoints, lastTrials[0]);
        lastTrialsPos[0] = insertInSorted(I[(size_t)lastTrials[0].nu - 1], lastTrials[0]);

        M = lastTrials[0].nu;
    }
    coincideX = false;

    double xNew, h;
    int t;
    vector<double> P(n), X(n);
    while (true) {
        // Steps 3, 4, 5, 6, Pre-7
        xNew = selectNewPoint(t);
        if (key == 3) {
            h = calcH(xNew);
            y(h, P);
        }

        // Step 7
        if (key != 3) {
            lastTrials[0] = newTrial(xNew);
        } else {
            x(P, hNu);
            if (!checkDensity(hNu[0])) {
                coincideX = true;
                numberFevals = this->numberFevals;
                return false;
            }
            lastTrials.clear();
            TrialConstrained trial = newTrial(hNu[0]);
            for (int i = 0; i < hNu.size(); i++) {
                trial.x = hNu[i];
                lastTrials.push_back(trial);
            }
        }

        numberTrials++;

        // Step 1
        for (int i = 0; i < lastTrials.size(); i++) {
            insertInSorted(trialPoints, lastTrials[i]);
        }

        // Step 2
        if (key != 3) {
            lastTrialsPos[0] = insertInSorted(I[(size_t)lastTrials[0].nu - 1], lastTrials[0]);
        } else {
            lastTrialsPos.resize(lastTrials.size());
            for (int i = 0; i < lastTrials.size(); i++) {
                lastTrialsPos[i] = insertInSorted(I[(size_t)lastTrials[0].nu - 1], lastTrials[i]);
            }
        }

        if (lastTrials[0].nu > M) {
            M = lastTrials[0].nu;
        }

        y(xNew, X);
        // Stop conditions
        if (euclideanDistance(X, XOpt) <= eps) {
            numberFevals = this->numberFevals;
            return true;
        }
        if (this->numberFevals >= maxFevals || numberTrials >= maxTrials) {
            numberFevals = this->numberFevals;
            return false;
        }
    }
}

bool MggsaMethod::solveTest(vector<double> XOpt, int &numberTrials, int &numberFevals) {
    return solveTest(XOpt, numberTrials, numberFevals, TypeSolve::SOLVE);
}
