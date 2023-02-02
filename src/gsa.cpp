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

inline int insert_in_sorted(vector<trial> &vec, trial tr) {
    vector<trial>::iterator iter = vec.begin();
    vector<trial>::iterator iterEnd = vec.end();
    int pos = 0;
    while(true) {
        if (iter == iterEnd || iter->x > tr.x) break;
        iter++; pos++;
    }
    vec.insert(iter, tr);
    return pos;
}

inline double search_min(vector<trial> &trials) {
    double z = numeric_limits<double>::infinity(), x;
    for (int i = 0; i < trials.size(); i++) {
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
    return (trial_points[t].x + trial_points[(size_t)t - 1].x) / 2 - (trial_points[t].z - trial_points[(size_t)t - 1].z) / (2 * m);
}

double gsa_method::selectNewPoint(int &t) {
    static double M = -1.0;

    // Step 2
    // with optimization(const)
    if (last_trial.x == B[0]) {
        M = abs((trial_points[1].z - trial_points[0].z) / (trial_points[1].x - trial_points[0].x));
    } else {
        M = max({M, abs((last_trial.z - trial_points[(size_t)last_trial_pos - 1].z) / 
                        (last_trial.x - trial_points[(size_t)last_trial_pos - 1].x)), 
                    abs((trial_points[(size_t)last_trial_pos + 1].z - last_trial.z) / 
                        (trial_points[(size_t)last_trial_pos + 1].x - last_trial.x))});
    }

    // without optimization
    // double M_tmp;
    // size_t size = trial_points.size();
    // M = 0.0;
    // for (int i = 1; i < size; i++) {
    //     M_tmp = abs(trial_points[i].z - trial_points[(size_t)i - 1].z) / (trial_points[i].x - trial_points[(size_t)i - 1].x);
    //     if (M_tmp > M) M = M_tmp;
    // }

    // Step 3
    m = (abs(M) <= epsilon) ? 1.0 : r * M;

    // Steps 4, 5
    double R = -numeric_limits<double>::infinity(), Rtmp = 0.0;
    double d_x;
    
    size_t size_tr_pt = trial_points.size();
    for (size_t i = 1; i < size_tr_pt; i++) {
        d_x = trial_points[i].x - trial_points[i - 1].x;
        Rtmp = m * d_x + pow(trial_points[i].z - trial_points[i - 1].z, 2) / (m * d_x) - 2 * (trial_points[i].z + trial_points[i - 1].z);
        if (Rtmp > R) {
            R = Rtmp;
            t = (int)i;
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
    last_trial = newTrial(B[0]);
    last_trial_pos = 1;
    while(true) {
        count++;

        // Steps 2, 3, 4, 5, 6
        x_k_1 = selectNewPoint(t);

        // Trial
        last_trial = newTrial(x_k_1);

        // Stop conditions
        if (trial_points[t].x - trial_points[(size_t)t - 1].x <= eps) {
            if (stop == Stop::ACCURACY || stop == Stop::ACCURNUMBER) break;
        }
        if (count >= Nmax) {
            if (stop == Stop::NUMBER || stop == Stop::ACCURNUMBER) break;
        }

        // Step 1
        last_trial_pos = insert_in_sorted(trial_points, last_trial);
    }
    x = search_min(trial_points);
}

void gsa_method::solve(int &count, vector<double> &X, Stop stop) {
    solve(count, X[0], stop);
}
