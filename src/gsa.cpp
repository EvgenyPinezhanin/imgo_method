#include <gsa.h>

#include <algorithm>
#if defined( _MSC_VER )
    #define _USE_MATH_DEFINES
    #include <math.h>
#else
    #include <cmath>
#endif

void addInSort(vector<trial> &vec, trial tr) {
    vector<trial>::iterator iter = vec.begin();
    vector<trial>::iterator iterEnd = vec.end();
    while(true) {
        if (iter == iterEnd || iter->x > tr.x) break;
        iter++;
    }
    vec.insert(iter, tr);
}

double searchMinX(vector<trial> &trials) {
    double z = trials[0].z;
    double x = trials[0].x;
    for (int i = 1; i < trials.size(); i++) {
        if (trials[i].z < z) {
            z = trials[i].z;
            x = trials[i].x;
        }
    }
    return x;
}

trial gsa_method::newTrial(double x) {
    return trial(x, f(x));
}

double gsa_method::newPoint(int t) {
    return (trial_points[t].x + trial_points[t - 1].x) / 2 - (trial_points[t].z - trial_points[t - 1].z) / (2 * m);
}

double gsa_method::selectNewPoint(int &t, trial last_trial) {
    static double M = -1.0;

    // Step 2
    if (last_trial.x == B[0]) {
        M = std::max({M, std::abs((last_trial.z - trial_points[t - 1].z) / (last_trial.x - trial_points[t - 1].x))});
    } else {
        M = std::max({M, std::abs((last_trial.z - trial_points[t - 1].z) / (last_trial.x - trial_points[t - 1].x)), 
                         std::abs((trial_points[t + 1].z - last_trial.z) / (trial_points[t + 1].x - last_trial.x))});
    }

    // Step 3
    m = (M == 0.0) ? 1.0 : r * M;

    // Steps 4, 5
    double d_x = trial_points[1].x - trial_points[0].x;
    double R = m * d_x + pow(trial_points[1].z - trial_points[0].z, 2) / (m * d_x) - 2 * (trial_points[1].z + trial_points[0].z);
    double Rtmp;
    t = 1;
    size_t size = trial_points.size();
    for (size_t i = 2; i < size; i++) {
        d_x = trial_points[i].x - trial_points[i - 1].x;
        Rtmp = m * d_x + pow(trial_points[i].z - trial_points[i - 1].z, 2) / (m * d_x) - 2 * (trial_points[i].z + trial_points[i - 1].z);
        if (Rtmp > R) {
            R = Rtmp;
            t = i;
        }
    }

    // Step 6
    return newPoint(t);
}

void gsa_method::solve(int &count, double &x, Stop stop) {
    trial_points.clear();

    trial_points.push_back(newTrial(A[0]));
    trial_points.push_back(newTrial(B[0]));
    count = 2;

    double x_k_1;
    int t = 1;
    trial tr = newTrial(B[0]);
    while(true) {
        x_k_1 = selectNewPoint(t, tr);
        tr = newTrial(x_k_1);

        // Step 1
        addInSort(trial_points, tr);

        count++;
        if (stop == ACCURACY) {
            if (trial_points[t].x - trial_points[t - 1].x <= eps) {
                break;
            }
        } else if (stop == NUMBER) {
            if (count >= Nmax) {
                break;
            }
        } else if (stop == ACCURNUMBER) {
            if (trial_points[t].x - trial_points[t - 1].x <= eps || count >= Nmax) {
                break;
            }
        }
    }
    x = searchMinX(trial_points);
}

void gsa_method::solve(int &count, vector<double> &X, Stop stop) {
    solve(count, X[0], stop);
}
