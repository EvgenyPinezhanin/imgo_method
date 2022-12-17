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

// #define DEBUG
// #define EPS
// #define TIME_TEST

#if defined(TIME_TEST)
    ofstream ofstr_test;
    vector<int> time_count;
    vector<double> time_mu, time_z_star, time_R, time_trial, time_add_trial, time_add_I;
    int start_time, end_time;
#endif

const double peano_a = 0.0, peano_b = 1.0, peano_random = 0.5;

inline double euclidean_distance(vector<double> val1, vector<double> val2) {
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

trial_constr mggsa_method::newTrial(double x) {
    trial_constr tr(x);
    vector<double> X(n);
    y(x, X);
    for (int j = 1; j <= m + 1; j++) {
        if ((f(X, j) > 0) || (j == m + 1)) {
            tr.z = f(X, j);
            tr.nu = j;
            break;
        }
    }

#if defined(DEBUG)
    cout << "trial: x = " << x << " X = " << X[0] << " Y = " << X[1] << " z = " << f(X, m + 1) << " nu = " << tr.nu << endl;
#endif

    return tr;
}

double mggsa_method::newPoint(int t) {
    if (trial_points[t].nu != trial_points[(size_t)t - 1].nu) {
        return (trial_points[t].x + trial_points[(size_t)t - 1].x) / 2.0;
    } else {
        return (trial_points[t].x + trial_points[(size_t)t - 1].x) / 2.0 -
                sgn(trial_points[t].z - trial_points[(size_t)t - 1].z) / (2.0 * r) *
                pow(abs(trial_points[t].z - trial_points[(size_t)t - 1].z) / mu[(size_t)trial_points[t].nu - 1], n);
    }
}

double mggsa_method::calc_h(double x_k_1) {
    double d = 1.0 / (pow(2.0, den * n) * (pow(2.0, n) - 1.0));
    double h = 0.0;
    while (x_k_1 - h > d) {
        h += d;
    }
    return h;
}

bool mggsa_method::check_density(double h_j) {
    for (int i = 0; i < trial_points.size(); i++) {
        if (abs(h_j - trial_points[i].x) <= 1e-14) return false; 
    }
    return true;
}

void mggsa_method::y(double x, vector<double> &X) {
    int d = (key != 3) ? den : den + 1;
    mapd(x, d, X.data(), n, key);
    for (int i = 0; i < n; i++) {
        X[i] = X[i] * (B[i] - A[i]) + (A[i] + B[i]) / 2.0;
    }
}

void mggsa_method::x(const vector<double> &P, vector<double> &X) {
    vector<double> P_correct(n);
    int size_x;
    X.resize(pow(2, n));
    for (int i = 0; i < n; i++) {
        P_correct[i] = (P[i] - (A[i] + B[i]) / 2.0) / (B[i] - A[i]);
    }
    invmad(den + 1, X.data(), pow(2, n), &size_x, P_correct.data(), n, incr);
    X.resize(size_x);
}

double mggsa_method::selectNewPoint(int &t) {
#if defined(TIME_TEST)
    start_time = clock();
#endif
    
    // Step 3
    // with optimization
    double mu_tmp;
    int nu_I = last_trials[0].nu - 1;
    size_t size_I = I[nu_I].size();
    for (int nu = 0; nu < m + 1; nu++) {
        if (!calc_I[nu]) mu[nu] = 0.0;
    }
    if (I[nu_I].size() >= 3) {
        for (int i = 0; i < last_trials_pos.size(); i++) {
            if (last_trials_pos[i] == 0) {
                mu[nu_I] = max({ mu[nu_I],
                                 abs(I[nu_I][last_trials_pos[i] + 1].z - I[nu_I][last_trials_pos[i]].z) / 
                                 pow(I[nu_I][last_trials_pos[i] + 1].x - I[nu_I][last_trials_pos[i]].x, 1.0 / n) });  
            } else if (last_trials_pos[i] == I[nu_I].size() - 1) {
                mu[nu_I] = max({ mu[nu_I],
                                 abs(I[nu_I][last_trials_pos[i]].z - I[nu_I][last_trials_pos[i] - 1].z) / 
                                 pow(I[nu_I][last_trials_pos[i]].x - I[nu_I][last_trials_pos[i] - 1].x, 1.0 / n) });
            } else {
                mu[nu_I] = max({ mu[nu_I],
                                 abs(I[nu_I][last_trials_pos[i]].z - I[nu_I][last_trials_pos[i] - 1].z) / 
                                 pow(I[nu_I][last_trials_pos[i]].x - I[nu_I][last_trials_pos[i] - 1].x, 1.0 / n),
                                 abs(I[nu_I][last_trials_pos[i] + 1].z - I[nu_I][last_trials_pos[i]].z) / 
                                 pow(I[nu_I][last_trials_pos[i] + 1].x - I[nu_I][last_trials_pos[i]].x, 1.0 / n) });
            }
        }
    } else if (I[nu_I].size() == 2) {
        mu[nu_I] = max({ mu[nu_I], abs(I[nu_I][1].z - I[nu_I][0].z) / pow(I[nu_I][1].x - I[nu_I][0].x, 1.0 / n) });
    }
    if (abs(mu[nu_I]) > 1e-14) calc_I[nu_I] = true;
    for (int nu = 0; nu < m + 1; nu++) {
        if (abs(mu[nu]) < 1e-14) mu[nu] = 1.0;
    }

    // without optimization
    // for (int nu = 0; nu < m + 1; nu++) {
    //     mu[nu] = 0.0;
    // }
    // double mu_tmp;
    // size_t size_I;
    // for (int nu = 0; nu < m + 1; nu++) {
    //     size_I = I[nu].size();
    //     for (int i = 1; i < size_I; i++) { // при i = 0 - нет j
    //         for (int j = 0; j < i; j++) {
    //             mu_tmp = abs(I[nu][i].z - I[nu][j].z) / pow(I[nu][i].x - I[nu][j].x, 1.0 / n);
    //             if (mu_tmp > mu[nu]) {
    //                 mu[nu] = mu_tmp;
    //             }
    //         }
    //     }
    //     if (abs(mu[nu]) < 1e-14) {
    //         mu[nu] = 1.0;
    //     };
    // }

#if defined(TIME_TEST)
    end_time = clock();
    time_mu.push_back((double)(end_time - start_time) / CLOCKS_PER_SEC);
#endif

#if defined(DEBUG)
    cout << "mu: ";
    for (int nu = 0; nu < m + 1; nu++) {
        cout << mu[nu] << " ";
    }
    cout << endl;
#endif

#if defined(TIME_TEST)
    start_time = clock();
#endif

    // Step 4
    for (int nu = 0; nu < m + 1; nu++) {
        if (I[nu].size() != 0) {
            if (nu + 1 != M) {
                z_star[nu] = -mu[nu] * d;
            } else {
                z_star[nu] = I[nu][0].z;
                size_I = I[nu].size();
                for (int i = 1; i < size_I; i++) {
                    if (I[nu][i].z < z_star[nu]) {
                        z_star[nu] = I[nu][i].z;
                    }
                }
            }
        }
    }

#if defined(TIME_TEST)
    end_time = clock();
    time_z_star.push_back((double)(end_time - start_time) / CLOCKS_PER_SEC);
#endif

#if defined(DEBUG)
    cout << "z_star: ";
    for (int nu = 0; nu < m + 1; nu++) {
        cout << z_star[nu] << " ";
    }
    cout << endl;
#endif

#if defined(TIME_TEST)
    start_time = clock();
#endif

    // Steps 5, 6
    double R = -numeric_limits<double>::infinity(), Rtmp = 0.0;
    double mu_v, z_star_v, d_x;

    size_t size_tr_pt = trial_points.size();
    for (size_t i = 1; i < size_tr_pt; i++) {
        d_x = pow(trial_points[i].x - trial_points[i - 1].x, 1.0 / n);
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

#if defined(DEBUG)
    cout << "R[" << trial_points[i - 1].x << ", " << trial_points[i].x << "] = " << Rtmp << endl;
#endif

    }

#if defined(TIME_TEST)
    end_time = clock();
    time_R.push_back((double)(end_time - start_time) / CLOCKS_PER_SEC);
#endif

    // Pre-Step 7
    return newPoint(t);
}

void mggsa_method::setM(int _m) {
    optimization_method_constrained::setM(_m);
    I.resize((size_t)m + 1);
    calc_I.resize((size_t)m + 1);
    mu.resize((size_t)m + 1);
    z_star.resize((size_t)m + 1);
}

void mggsa_method::getPoints(vector<vector<double>> &points_vec) {
    points_vec.clear();
    vector<double> point(n);
    for (int i = 0; i < trial_points.size(); i++) {
        y(trial_points[i].x, point);
        points_vec.push_back(point);
    }
}

void mggsa_method::solve(int &count, vector<double> &X, Stop stop) {
#if defined(TIME_TEST)
    ofstr_test.open("output_data/mggsa_time_test.txt");
    if (!ofstr_test.is_open()) cerr << "File opening error\n";
    time_count.clear();
    time_mu.clear();
    time_z_star.clear();
    time_R.clear();
    time_trial.clear();
    time_add_trial.clear();
    time_add_I.clear();
#endif

    for (int nu = 0; nu < m + 1; nu++) {
        I[nu].clear();
        calc_I[nu] = false;
    }
    trial_points.clear();
    X.resize(n);

    last_trials[0] = trial_constr{peano_a, -1.0, 0};
    trial_points.push_back(last_trials[0]);
    last_trials[0].x = peano_b;
    trial_points.push_back(last_trials[0]);

#if defined(DEBUG)
    cout << "Trial number: " << 1 << endl;
#endif

#if defined(TIME_TEST)
    start_time = clock();
#endif

    last_trials[0] = newTrial(peano_random);

#if defined(TIME_TEST)
    end_time = clock();
#endif

#if defined(TIME_TEST)
    time_count.push_back(0);
    time_mu.push_back(0.0);
    time_z_star.push_back(0.0);
    time_R.push_back(0.0);
    time_trial.push_back((double)(end_time - start_time) / CLOCKS_PER_SEC);
#endif

#if defined(TIME_TEST)
    start_time = clock();
#endif

    insert_in_sorted(trial_points, last_trials[0]);

#if defined(TIME_TEST)
    end_time = clock();
    time_add_trial.push_back((double)(end_time - start_time) / CLOCKS_PER_SEC);
#endif

#if defined(TIME_TEST)
    start_time = clock();
#endif

    last_trials_pos[0] = insert_in_sorted(I[(size_t)last_trials[0].nu - 1], last_trials[0]);

#if defined(TIME_TEST)
    end_time = clock();
    time_add_I.push_back((double)(end_time - start_time) / CLOCKS_PER_SEC);
#endif

    M = last_trials[0].nu;
    count = 1;

    double x_k_1, h, delta_t;
    int t = 1;
    vector<double> P(n);
    while(true) {
    #if defined(TIME_TEST)
        time_count.push_back(count);
    #endif

    #if defined(DEBUG)
        cout << "Trial number: " << count + 1 << endl;
    #endif

        // Steps 3, 4, 5, 6, Pre-7
        x_k_1 = selectNewPoint(t);
        if (key == 3) {
            h = calc_h(x_k_1);
            y(h, P);
        }

        delta_t = pow(trial_points[t].x - trial_points[(size_t)t - 1].x, 1.0 / n);

    #if defined(TIME_TEST)
        start_time = clock();
    #endif

        // Step 7
        if (key != 3) {
            last_trials.resize(1);
            last_trials[0] = newTrial(x_k_1);
        } else {
            x(P, h_nu);
            if (!check_density(h_nu[0])) break;
            last_trials.clear();
            for (int i = 0; i < h_nu.size(); i++) {
                last_trials.push_back(newTrial(h_nu[i]));
            }
        }

    #if defined(TIME_TEST)
        end_time = clock();
        time_trial.push_back((double)(end_time - start_time) / CLOCKS_PER_SEC);
    #endif
        
    #if defined(TIME_TEST)
        start_time = clock();
    #endif

        // Step 1
        for (int i = 0; i < last_trials.size(); i++) {
            insert_in_sorted(trial_points, last_trials[i]);
        }

    #if defined(TIME_TEST)
        end_time = clock();
        time_add_trial.push_back((double)(end_time - start_time) / CLOCKS_PER_SEC);
    #endif

    #if defined(TIME_TEST)
        start_time = clock();
    #endif

        // Step 2
        if (key != 3) {
            last_trials_pos.resize(1);
            last_trials_pos[0] = insert_in_sorted(I[(size_t)last_trials[0].nu - 1], last_trials[0]);
        } else {
            last_trials_pos.resize(last_trials.size());
            for (int i = 0; i < last_trials.size(); i++) {
                last_trials_pos[i] = insert_in_sorted(I[(size_t)last_trials[0].nu - 1], last_trials[i]);
            }
        }
        
    #if defined(TIME_TEST)
        end_time = clock();
        time_add_I.push_back((double)(end_time - start_time) / CLOCKS_PER_SEC);
    #endif

        if (last_trials[0].nu > M) {
            M = last_trials[0].nu;
        }

        count++;

        #if defined(EPS)
            cout << pow(abs(trial_points[t].x - trial_points[t - 1].x), 1.0 / n) << endl;
        #endif

        if (delta_t <= eps) {
            if (stop == Stop::ACCURACY || stop == Stop::ACCURNUMBER) break;
        }
        if (count >= Nmax) {
            if (stop == Stop::NUMBER || stop == Stop::ACCURNUMBER) break;
        }
    }
    y(search_min(trial_points, m), X);

#if defined(TIME_TEST)
    for (int i = 0; i < count - 1; i++) {
        ofstr_test << time_count[i] << " " << time_mu[i] << " " << time_z_star[i] << " "
                   << time_R[i] << " " << time_trial[i] << " "  << time_add_trial[i] << " " << time_add_I[i] << endl;
    }
    ofstr_test.close();

    // Plotting the functions of time(works with gnuplot)
    int error;
    setenv("QT_QPA_PLATFORM", "xcb", false);
    error = system("chmod +x scripts/time_test.gp");
    if (error != 0) {
        cerr << "Error chmod" << endl;
    }
    error = system("gnuplot -c scripts/time_test.gp");
    if (error != 0) {
        cerr << "Error gnuplot" << endl;
    }
#endif
}

bool mggsa_method::solve_test(vector<double> X_opt, int &count, Stop stop) {
    for (int nu = 0; nu < m + 1; nu++) {
        I[nu].clear();
        calc_I[nu] = false;
    }
    trial_points.clear();

    last_trials[0] = trial_constr{peano_a, -1.0, 0};
    trial_points.push_back(last_trials[0]);
    last_trials[0].x = peano_b;
    trial_points.push_back(last_trials[0]);

    last_trials[0] = newTrial(peano_random);
    insert_in_sorted(trial_points, last_trials[0]);
    last_trials_pos[0] = insert_in_sorted(I[(size_t)last_trials[0].nu - 1], last_trials[0]);

    M = last_trials[0].nu;
    count = 1;

    double x_k_1, h, delta_t;
    int t = 1;
    vector<double> P(n), X(n);
    while (true) {
        // Steps 3, 4, 5, 6, Pre-7
        x_k_1 = selectNewPoint(t);
        if (key == 3) {
            h = calc_h(x_k_1);
            y(h, P);
        }

        // Step 7
        if (key != 3) {
            last_trials.resize(1);
            last_trials[0] = newTrial(x_k_1);
        } else {
            x(P, h_nu);
            if (!check_density(h_nu[0])) return false;
            last_trials.clear();
            for (int i = 0; i < h_nu.size(); i++) {
                last_trials.push_back(newTrial(h_nu[i]));
            }
        }

        // Step 1
        for (int i = 0; i < last_trials.size(); i++) {
            insert_in_sorted(trial_points, last_trials[i]);
        }

        // Step 2
        if (key != 3) {
            last_trials_pos.resize(1);
            last_trials_pos[0] = insert_in_sorted(I[(size_t)last_trials[0].nu - 1], last_trials[0]);
        } else {
            last_trials_pos.resize(last_trials.size());
            for (int i = 0; i < last_trials.size(); i++) {
                last_trials_pos[i] = insert_in_sorted(I[(size_t)last_trials[0].nu - 1], last_trials[i]);
            }
        }

        if (last_trials[0].nu > M) {
            M = last_trials[0].nu;
        }

        count++;

        y(x_k_1, X);
        if (euclidean_distance(X, X_opt) <= eps) {
            if (stop == Stop::ACCURACY || stop == Stop::ACCURNUMBER) return true;
        }
        if (count >= Nmax) {
            if (stop == Stop::NUMBER || stop == Stop::ACCURNUMBER) return false;
        }
    }
}
