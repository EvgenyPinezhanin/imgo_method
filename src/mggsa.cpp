#include <mggsa.h>

#include <iostream>
#include <fstream>
#include <ctime>
#if defined( _MSC_VER )
    #define _USE_MATH_DEFINES
    #include <math.h>
#else
    #include <cmath>
#endif

#include <peano.h>

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

double euclidean_distance(vector<double> val1, vector<double> val2) {
    double res = 0.0;
    size_t size = val1.size();
    for (int i = 0; i < size; i++) {
        res += (val1[i] - val2[i]) * (val1[i] - val2[i]);
    }
    res = sqrt(res);
    return res;
}

template<typename T> 
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

void addInSort(vector<trial_constr> &vec, trial_constr tr) {
    vector<trial_constr>::iterator iter = vec.begin();
    vector<trial_constr>::iterator iterEnd = vec.end();
    while(true) {
        if (iter == iterEnd || iter->x > tr.x) break;
        iter++;
    }
    vec.insert(iter, tr);
}

double searchMinXTrial(vector<trial_constr> &trials, int m) {
    double z = 0.0, x = 0.0;
    int k = 0;
    for (int i = 0; i < trials.size(); i++) {
        if (trials[i].nu == m + 1) {
            z = trials[i].z;
            x = trials[i].x;
            k = i + 1;
            break;
        }
    }
    for (int i = k; i < trials.size(); i++) {
        if (trials[i].z < z && trials[i].nu == m + 1) {
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
    cout << x << " " << X[0] << " " << X[1] << " " << f(X, m + 1) << " " << tr.nu << endl;
#endif

    return tr;
}

double mggsa_method::newPoint(int t) {
    if (trial_points[t].nu != trial_points[t - 1].nu) {
        return (trial_points[t].x + trial_points[t - 1].x) / 2.0;
    } else {
        return (trial_points[t].x + trial_points[t - 1].x) / 2.0 - 
                sgn(trial_points[t].z - trial_points[t - 1].z) / (2.0 * r) * 
                pow(abs(trial_points[t].z - trial_points[t - 1].z) / mu[trial_points[t].nu - 1], n);
    }
}

double mggsa_method::selectNewPoint(int &t, trial_constr last_trial) {
#if defined(TIME_TEST)
    start_time = clock();
#endif
    
    // Step 3
    double mu_tmp;
    int nu_I = last_trial.nu - 1;
    size_t size_I = I[nu_I].size();
    for (int nu = 0; nu < m + 1; nu++) {
        if (!calc_I[nu]) mu[nu] = 0.0;
    }
    for (int i = 0; i < size_I; i++) {
        if (I[last_trial.nu - 1][i].x != last_trial.x) {
            mu_tmp = abs(I[nu_I][i].z - last_trial.z) / pow(abs(I[last_trial.nu - 1][i].x - last_trial.x), 1.0 / n);
            if (mu_tmp > mu[nu_I]) {
                mu[nu_I] = mu_tmp;
                if (abs(mu[nu_I]) > 1e-14) calc_I[nu_I] = true;
            }
        }
    }
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
    //     if (abs(mu[nu]) < 6e-13) {
    //         mu[nu] = 1.0;
    //     };
    // }

    // not work
    // mu[last_trial.nu - 1] = 0.0;
    // for (int i = 1; i < size_I; i++) {
    //     mu_tmp = abs(I[last_trial.nu - 1][i].z - I[last_trial.nu - 1][i - 1].z) / deltaX(I[last_trial.nu - 1][i - 1].x, I[last_trial.nu - 1][i].x);
    //     if (mu_tmp > mu[last_trial.nu - 1]) {
    //         mu[last_trial.nu - 1] = mu_tmp;
    //     }
    // }
    // if (abs(mu[last_trial.nu - 1]) < 6e-13) {
    //     mu[last_trial.nu - 1] = 1.0;
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
    t = 1;
    double d_x = pow(trial_points[1].x - trial_points[0].x, 1.0 / n); 
    double R, Rtmp = 0.0;
    double mu_v, z_star_v;
    if (trial_points[1].nu == trial_points[0].nu) {
        mu_v = mu[trial_points[1].nu - 1];
        z_star_v = z_star[trial_points[1].nu - 1];
        R = d_x + pow(trial_points[1].z - trial_points[0].z, 2) / (r * r * mu_v * mu_v * d_x) -
            2.0 * (trial_points[1].z + trial_points[0].z - 2.0 * z_star_v) / (r * mu_v);
    } else if (trial_points[0].nu < trial_points[1].nu) {
        mu_v = mu[trial_points[1].nu - 1];
        z_star_v = z_star[trial_points[1].nu - 1];
        R = 2.0 * d_x  - 4.0 * (trial_points[1].z - z_star_v) / (r * mu_v);
    } else {
        mu_v = mu[trial_points[0].nu - 1];
        z_star_v = z_star[trial_points[0].nu - 1];
        R = 2.0 * d_x  - 4.0 * (trial_points[0].z - z_star_v) / (r * mu_v);
    }

#if defined(DEBUG)
    cout << "R[" << trial_points[0].x << ", " << trial_points[1].x << "] = " << R << endl;
#endif

    size_t size_tr_pt = trial_points.size();
    for (size_t i = 2; i < size_tr_pt; i++) {
        d_x = pow(trial_points[i].x - trial_points[i - 1].x, 1.0 / n);
        if (trial_points[i].nu == trial_points[i - 1].nu) {
            mu_v = mu[trial_points[i].nu - 1];
            z_star_v = z_star[trial_points[i].nu - 1];
            Rtmp = d_x + pow(trial_points[i].z - trial_points[i - 1].z, 2) / (r * r * mu_v * mu_v * d_x) -
                   2.0 * (trial_points[i].z + trial_points[i - 1].z - 2.0 * z_star_v) / (r * mu_v);
        } else if (trial_points[i - 1].nu < trial_points[i].nu) {
            mu_v = mu[trial_points[i].nu - 1];
            z_star_v = z_star[trial_points[i].nu - 1];
            Rtmp = 2.0 * d_x  - 4.0 * (trial_points[i].z - z_star_v) / (r * mu_v);
        } else  {
            mu_v = mu[trial_points[i - 1].nu - 1];
            z_star_v = z_star[trial_points[i - 1].nu - 1];
            Rtmp = 2.0 * d_x  - 4.0 * (trial_points[i - 1].z - z_star_v) / (r * mu_v);
        }
        if (Rtmp > R) {
            R = Rtmp;
            t = i;
        }

#if defined(DEBUG)
    cout << "R[" << trial_points[i - 1].x << ", " << trial_points[i].x << "] = " << Rtmp << endl;
#endif

    }

#if defined(TIME_TEST)
    end_time = clock();
    time_R.push_back((double)(end_time - start_time) / CLOCKS_PER_SEC);
#endif

    // Step 7
    return newPoint(t);
}

void mggsa_method::setM(int _m) {
    optimization_method_constrained::setM(_m);
    I.resize(m + 1);
    calc_I.resize(m + 1);
    mu.resize(m + 1);
    z_star.resize(m + 1);
}

void mggsa_method::getPoints(vector<vector<double>> &points_vec) {
    points_vec.clear();
    vector<double> point(n);
    for (int i = 0; i < trial_points.size(); i++) {
        y(trial_points[i].x, point);
        points_vec.push_back(point);
    }
}

void mggsa_method::y(double x, vector<double> &X) {
    mapd(x, den, X.data(), n, key);
    for (int i = 0; i < n; i++) {
        X[i] = X[i] * (B[i] - A[i]) + (A[i] + B[i]) / 2.0;
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
    M = 0;

    trial_constr tr{peano_a, -1.0, 0};
    trial_points.push_back(tr);
    tr.x = peano_b;
    trial_points.push_back(tr);

#if defined(TIME_TEST)
    start_time = clock();
#endif

    tr = newTrial(peano_random);

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

    addInSort(trial_points, tr);

#if defined(TIME_TEST)
    end_time = clock();
    time_add_trial.push_back((double)(end_time - start_time) / CLOCKS_PER_SEC);
#endif

#if defined(TIME_TEST)
    start_time = clock();
#endif

    addInSort(I[tr.nu - 1], tr);

#if defined(TIME_TEST)
    end_time = clock();
    time_add_I.push_back((double)(end_time - start_time) / CLOCKS_PER_SEC);
#endif

    if (tr.nu > M) {
        M = tr.nu;
    }
    count = 1;

    double x_k_1;
    int t = 1;
    while(true) {
    #if defined(TIME_TEST)
        time_count.push_back(count);
    #endif

        x_k_1 = selectNewPoint(t, tr);

    #if defined(TIME_TEST)
        start_time = clock();
    #endif

        tr = newTrial(x_k_1);

    #if defined(TIME_TEST)
        end_time = clock();
        time_trial.push_back((double)(end_time - start_time) / CLOCKS_PER_SEC);
    #endif
        
    #if defined(TIME_TEST)
        start_time = clock();
    #endif

        // Step 1
        addInSort(trial_points, tr);

    #if defined(TIME_TEST)
        end_time = clock();
        time_add_trial.push_back((double)(end_time - start_time) / CLOCKS_PER_SEC);
    #endif

    #if defined(TIME_TEST)
        start_time = clock();
    #endif

        // Step 2
        addInSort(I[tr.nu - 1], tr);

    #if defined(TIME_TEST)
        end_time = clock();
        time_add_I.push_back((double)(end_time - start_time) / CLOCKS_PER_SEC);
    #endif

        if (tr.nu > M) {
            M = tr.nu;
        }

        count++;

        #if defined(EPS)
            cout << pow(abs(trial_points[t].x - trial_points[t - 1].x), 1.0 / n) << endl;
        #endif

        if (pow(trial_points[t].x - trial_points[t - 1].x, 1.0 / n) <= eps) {
            if (stop == ACCURACY || stop == ACCURNUMBER) break;
        }
        if (count >= Nmax) {
            if (stop == NUMBER || stop == ACCURNUMBER) break;
        }
    }
    y(searchMinXTrial(trial_points, m), X);

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
    error = system("gnuplot -p -c scripts/time_test.gp");
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
    Nmax = count;
    M = 0;

    trial_constr tr{peano_a, -1.0, 0};
    trial_points.push_back(tr);
    tr.x = peano_b;
    trial_points.push_back(tr);
    tr = newTrial(peano_random);
    addInSort(trial_points, tr);
    addInSort(I[tr.nu - 1], tr);
    if (tr.nu > M) {
        M = tr.nu;
    }
    count = 1;

    double x_k_1;
    int t = 1;
    vector<double> X(n);
    while (true) {
        x_k_1 = selectNewPoint(t, tr);
        tr = newTrial(x_k_1);

        // Step 1
        addInSort(trial_points, tr);

        // Step 2
        addInSort(I[tr.nu - 1], tr);

        if (tr.nu > M) {
            M = tr.nu;
        }

        count++;
        y(x_k_1, X);
        if (euclidean_distance(X, X_opt) <= eps) {
            if (stop == ACCURACY || stop == ACCURNUMBER) return true;
        }
        if (count >= Nmax) {
            if (stop == NUMBER || stop == ACCURNUMBER) return false;
        }
    }
}
