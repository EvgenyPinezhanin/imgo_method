#include <imgo.h>

#include <algorithm>
#if defined( _MSC_VER )
    #define _USE_MATH_DEFINES
    #include <math.h>
#else
    #include <cmath>
#endif

inline int insert_in_sorted(vector<trial_constr> &vec, trial_constr tr) {
    vector<trial_constr>::iterator iter = vec.begin();
    vector<trial_constr>::iterator iterEnd = vec.end();
    int pos = 0;
    while(true) {
        if (iter == iterEnd || iter->x > tr.x) break;
        iter++; pos++;
    }
    vec.insert(iter, tr);
    return pos;
}

inline double search_min(vector<trial_constr> &trials, int m) {
    double z = numeric_limits<double>::infinity(), x;
    for (int i = 0; i < trials.size(); i++) {
        if (trials[i].nu == m + 1 && trials[i].z < z) {
            z = trials[i].z;
            x = trials[i].x;
        }
    }
    return x;
}

trial_constr imgo_method::newTrial(double x) {
    trial_constr tr(x);
    for (int j = 1; j <= m + 1; j++) {
        if ((f(x, j) > 0) || (j == m + 1)) {
            tr.z = f(x, j);
            tr.nu = j;
            break;
        }
    }
    return tr;
}

double imgo_method::newPoint(int t) {
    if (trial_points[t].nu != trial_points[(size_t)t - 1].nu) {
        return (trial_points[t].x + trial_points[(size_t)t - 1].x) / 2.0;
    } else {
        return (trial_points[t].x + trial_points[(size_t)t - 1].x) / 2.0 - 
               (trial_points[t].z - trial_points[(size_t)t - 1].z) / (2.0 * r * mu[(size_t)trial_points[t].nu - 1]);
    }
}

double imgo_method::selectNewPoint(int &t) {
    // Step 3
    // with optimization
    double mu_tmp;
    int nu_I = last_trial.nu - 1;
    size_t size_I = I[nu_I].size();
    for (int nu = 0; nu < m + 1; nu++) {
        if (!calc_I[nu]) mu[nu] = 0.0;
    }
    if (I[nu_I].size() >= 3) {
        if (last_trial_pos == 0) {
            mu[nu_I] = max({ mu[nu_I],
                             abs(I[nu_I][last_trial_pos + 1].z - I[nu_I][last_trial_pos].z) / 
                             pow(I[nu_I][last_trial_pos + 1].x - I[nu_I][last_trial_pos].x, 1.0 / n) });  
        } else if (last_trial_pos == I[nu_I].size() - 1) {
            mu[nu_I] = max({ mu[nu_I],
                             abs(I[nu_I][last_trial_pos].z - I[nu_I][last_trial_pos - 1].z) / 
                             pow(I[nu_I][last_trial_pos].x - I[nu_I][last_trial_pos - 1].x, 1.0 / n) });
        } else {
            mu[nu_I] = max({ mu[nu_I],
                             abs(I[nu_I][last_trial_pos].z - I[nu_I][last_trial_pos - 1].z) / 
                             pow(I[nu_I][last_trial_pos].x - I[nu_I][last_trial_pos - 1].x, 1.0 / n),
                             abs(I[nu_I][last_trial_pos + 1].z - I[nu_I][last_trial_pos].z) / 
                             pow(I[nu_I][last_trial_pos + 1].x - I[nu_I][last_trial_pos].x, 1.0 / n) });
        }
    } else if (I[nu_I].size() == 2) {
        mu[nu_I] = max({ mu[nu_I], abs(I[nu_I][1].z - I[nu_I][0].z) / pow(I[nu_I][1].x - I[nu_I][0].x, 1.0 / n) });
    }
    if (abs(mu[nu_I]) > 1e-14) calc_I[nu_I] = true;
    for (int nu = 0; nu < m + 1; nu++) {
        if (abs(mu[nu]) < 1e-14) mu[nu] = 1.0;
    }

    // without optimization
    // double mu_tmp;
    // int nu_I = last_trial.nu - 1;
    // size_t size_I = I[nu_I].size();
    // for (int nu = 0; nu < m + 1; nu++) {
    //     if (!calc_I[nu]) mu[nu] = 0.0;
    // }
    // for (int i = 1; i < size_I; i++) {
    //     mu_tmp = abs(I[nu_I][i].z - I[nu_I][(size_t)i - 1].z) / (I[nu_I][i].x - I[nu_I][(size_t)i - 1].x);
    //     if (mu_tmp > mu[nu_I]) {
    //         mu[nu_I] = mu_tmp;
    //         if (abs(mu[nu_I]) > 1e-14) calc_I[nu_I] = true;
    //     }
    // }
    // for (int nu = 0; nu < m + 1; nu++) {
    //     if (abs(mu[nu]) < 1e-14) mu[nu] = 1.0;
    // }

    // Step 4
    for (int nu = 0; nu < m + 1; nu++) {
        if (I[nu].size() != 0) {
            z_star[nu] = I[nu][0].z;
            size_I = I[nu].size();
            for (int i = 1; i < size_I; i++) {
                if (I[nu][i].z < z_star[nu]) {
                    z_star[nu] = I[nu][i].z;
                }
            }
            if (z_star[nu] <= 0.0 && nu != m) {
                z_star[nu] = -mu[nu] * d;
            }
        }
    }

    // Steps 5, 6
    double R = -numeric_limits<double>::infinity(), Rtmp = 0.0;
    double mu_v, z_star_v, d_x;

    size_t size_tr_pt = trial_points.size();
    for (size_t i = 1; i < size_tr_pt; i++) {
        d_x = trial_points[i].x - trial_points[i - 1].x;
        if (trial_points[i].nu == trial_points[i - 1].nu) {
            mu_v = mu[(size_t)trial_points[i].nu - 1];
            z_star_v = z_star[(size_t)trial_points[i].nu - 1];
            Rtmp = d_x + pow(trial_points[i].z - trial_points[i - 1].z, 2) / (r * r * mu_v * mu_v * d_x) -
                   2.0 * (trial_points[i].z + trial_points[i - 1].z - 2.0 * z_star_v) / (r * mu_v);
        } else if (trial_points[i - 1].nu < trial_points[i].nu) {
            mu_v = mu[(size_t)trial_points[i].nu - 1];
            z_star_v = z_star[(size_t)trial_points[i].nu - 1];
            Rtmp = 2.0 * d_x  - 4.0 * (trial_points[i].z - z_star_v) / (r * mu_v);
        } else  {
            mu_v = mu[(size_t)trial_points[i - 1].nu - 1];
            z_star_v = z_star[(size_t)trial_points[i - 1].nu - 1];
            Rtmp = 2.0 * d_x  - 4.0 * (trial_points[i - 1].z - z_star_v) / (r * mu_v);
        }
        if (Rtmp > R) {
            R = Rtmp;
            t = (int)i;
        }
    }

    // Step 7
    return newPoint(t);
}

void imgo_method::setM(int _m) {
    optimization_method_constrained::setM(_m);
    I.resize((size_t)m + 1);
    calc_I.resize((size_t)m + 1);
    mu.resize((size_t)m + 1);
    z_star.resize((size_t)m + 1);
}

void imgo_method::solve(int &count, double &x, Stop stop) {
    for (int i = 0; i < I.size(); i++) {
        I[i].clear();
        calc_I[i] = false;
    }
    trial_points.clear();

    last_trial = newTrial(A[0]);
    trial_points.push_back(last_trial);
    insert_in_sorted(I[(size_t)last_trial.nu - 1], last_trial);
    last_trial = newTrial(B[0]);
    trial_points.push_back(last_trial);
    last_trial_pos = insert_in_sorted(I[(size_t)last_trial.nu - 1], last_trial);
    count = 2;

    double x_k_1;
    int t = 1;
    while(true) {
        x_k_1 = selectNewPoint(t);
        last_trial = newTrial(x_k_1);

        // Step 1
        insert_in_sorted(trial_points, last_trial);

        // Step 2
        last_trial_pos = insert_in_sorted(I[(size_t)last_trial.nu - 1], last_trial);

        count++;
        if (trial_points[t].x - trial_points[(size_t)t - 1].x <= eps) {
            if (stop == Stop::ACCURACY || stop == Stop::ACCURNUMBER) break;
        }
        if (count >= Nmax) {
            if (stop == Stop::NUMBER || stop == Stop::ACCURNUMBER) break;
        }
    }
    x = search_min(trial_points, m);
}

void imgo_method::solve(int &count, vector<double> &X, Stop stop) {
    solve(count, X[0], stop);
}

bool imgo_method::solve_test(double x_opt, int &count, Stop stop) {
    for (int nu = 0; nu < m + 1; nu++) {
        I[nu].clear();
        calc_I[nu] = false;
    }
    trial_points.clear();
    Nmax = count;

    last_trial = newTrial(A[0]);
    trial_points.push_back(last_trial);
    insert_in_sorted(I[(size_t)last_trial.nu - 1], last_trial);
    last_trial = newTrial(B[0]);
    trial_points.push_back(last_trial);
    last_trial_pos = insert_in_sorted(I[(size_t)last_trial.nu - 1], last_trial);
    count = 2;

    double x_k_1;
    int t = 1;
    while (true) {
        x_k_1 = selectNewPoint(t);
        last_trial = newTrial(x_k_1);

        // Step 1
        insert_in_sorted(trial_points, last_trial);

        // Step 2
        last_trial_pos = insert_in_sorted(I[(size_t)last_trial.nu - 1], last_trial);

        count++;
        if (abs(x_k_1 - x_opt) <= eps) {
            if (stop == Stop::ACCURACY || stop == Stop::ACCURNUMBER) return true;
        }
        if (count >= Nmax) {
            if (stop == Stop::NUMBER || stop == Stop::ACCURNUMBER) return false;
        }
    }
}
