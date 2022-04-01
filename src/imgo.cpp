#include<imgo.h>
#include<iostream>

//#define DEBUG

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

imgo_method::imgo_method(double (*_f)(double, int), size_t _m, double _a, double _b, double _eps, double _r, double _d) 
    : f(_f), m(_m), a(_a), b(_b), eps(_eps), r(_r), d(_d), n(1), Nmax(0) {
    I.resize(m + 1);
    mu.resize(m + 1);
    z_star.resize(m + 1);
}

imgo_method::imgo_method(double (*_f)(vector<double>, int), int _n, size_t _m, vector<double> _A, 
                         vector<double> _B, double _eps, double _r, double _d, int _den, int _key)
    : f_md(_f), n(_n), m(_m), A(_A), B(_B), eps(_eps), r(_r), d(_d), den(_den), key(_key), M(0), Nmax(0) {
    I.resize(m + 1);
    mu.resize(m + 1);
    z_star.resize(m + 1);
}


void imgo_method::addInSort(vector<trial> &vec, trial tr) {
    vector<trial>::iterator iter = vec.begin();
    vector<trial>::iterator iterEnd = vec.end();
    while(true) {
        if (iter == iterEnd || iter->x > tr.x) break;
        iter++;
    }
    vec.insert(iter, tr);
}

double imgo_method::searchMinX() {
    double z = 0.0;
    double x = 0.0;
    int k = 0;
    for (int i = 0; i < trial_points.size(); i++) {
        if (trial_points[i].nu == m + 1) {
            z = trial_points[i].z;
            x = trial_points[i].x;
            k = i + 1;
            break;
        }
    }
    for (int i = k; i < trial_points.size(); i++) {
        if (trial_points[i].z < z && trial_points[i].nu == m + 1) {
            z = trial_points[i].z;
            x = trial_points[i].x;
        }
    }
    return x;
}

trial imgo_method::newTrial(double x) {
    int nu = 0;
    double z = 0.0;
    for (int j = 1; j <= m + 1; j++) {
        if ((f(x, j) > 0) || (j == m + 1)) {
            z = f(x, j);
            nu = j;
            break;
        }
    }
    trial tr{x, z, nu};
    return tr;
}

trial imgo_method::newTrial_md(double x) {
    vector<double> X(n);
    int nu = 0;
    double z = 0.0;
    y(x, X);
    for (int j = 1; j <= m + 1; j++) {
        if ((f_md(X, j) > 0) || (j == m + 1)) {
            z = f_md(X, j);
            nu = j;
            break;
        }
    }

#if defined(DEBUG)
    std::cout << x << " " << X[0] << " " << X[1] << " " << f_md(X, m + 1) << " " << nu << std::endl;
#endif

    trial tr{x, z, nu};
    return tr;
}

double imgo_method::deltaX(double x1, double x2) {
    if (n == 1) {
        return x2 - x1;
    } else {
        return pow(x2 - x1, 1.0 / n);
    }
}

double imgo_method::newPoint(size_t t) {
    if (trial_points[t].nu != trial_points[t - 1].nu) {
        return (trial_points[t].x + trial_points[t - 1].x) / 2.0;
    } else {
        if (n == 1) {
            return (trial_points[t].x + trial_points[t - 1].x) / 2.0 
                   - (trial_points[t].z - trial_points[t - 1].z) / (2.0 * r * mu[trial_points[t].nu - 1]);
        } else {
            return (trial_points[t].x + trial_points[t - 1].x) / 2.0 
                   - sgn(trial_points[t].z - trial_points[t - 1].z) / (2.0 * r) 
                   * pow(abs(trial_points[t].z - trial_points[t - 1].z) / mu[trial_points[t].nu - 1], n);
        }

    }
}

double imgo_method::selectNewPoint(size_t &t) {
    // Шаг 3
    for (int nu = 0; nu < m + 1; nu++) {
        mu[nu] = 0.0;
    }
    double mu_tmp;
    size_t size_I;
    for (int nu = 0; nu < m + 1; nu++) {
        size_I = I[nu].size();
        for (int i = 1; i < size_I; i++) { // при i = 0 - нет j
            for (int j = 0; j < i; j++){
                mu_tmp = abs(I[nu][i].z - I[nu][j].z) / deltaX(I[nu][j].x, I[nu][i].x);
                if (mu_tmp > mu[nu]) {
                    mu[nu] = mu_tmp;
                }
            }
        }
        if (abs(mu[nu]) < 6e-13) {
            mu[nu] = 1.0;
        };
    }

#if defined(DEBUG)
    cout << "mu: ";
    for (int nu = 0; nu < m + 1; nu++) {
        cout << mu[nu] << " ";
    }
    cout << endl;
#endif

    // Шаг 4
    if (n >= 2) {
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
    } else {
        for (int nu = 0; nu < m + 1; nu++) {
            if (I[nu].size() != 0) {
                z_star[nu] = I[nu][0].z;
                size_I = I[nu].size();
                for (int i = 1; i < size_I; i++) {
                    if (I[nu][i].z < z_star[nu]) {
                        z_star[nu] = I[nu][i].z;
                    }
                }
                if (z_star[nu] <= 0 && nu != m) {
                    z_star[nu] = -mu[nu] * d;
                }
            }
        }
    }

#if defined(DEBUG)
    cout << "z_star: ";
    for (int nu = 0; nu < m + 1; nu++) {
        cout << z_star[nu] << " ";
    }
    cout << endl;
#endif 

    // Шаг 5, 6
    t = 1;
    double d_x = deltaX(trial_points[0].x, trial_points[1].x); 
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
        d_x = deltaX(trial_points[i - 1].x, trial_points[i].x);
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

    // Шаг 7
    return newPoint(t);
}

void imgo_method::setFunc(double (*_f)(double, int)) {
    f = _f;
}

void imgo_method::setFunc(double (*_f_md)(vector<double>, int)) {
    f_md = _f_md;
}

void imgo_method::setA(double _a) {
    a = _a;
}

void imgo_method::setA(vector<double> _A) {
    A = _A;
}

void imgo_method::setB(double _b) {
    b = _b;
}

void imgo_method::setB(vector<double> _B) {
    B = _B;
}

void imgo_method::setM(size_t _m) {
    m = _m;
    I.resize(m + 1);
    mu.resize(m + 1);
    z_star.resize(m + 1);
}

void imgo_method::setEps(double _eps) {
    eps = _eps;
}

void imgo_method::setR(double _r) {
    r = _r;
}

void imgo_method::setD(double _d) {
    d = _d;
}

void imgo_method::setN(int _n) {
    n = _n;
}

void imgo_method::setDen(int _den) {
    den = _den;
}

void imgo_method::setKey(int _key) {
    key = _key;
}

void imgo_method::getTrialPoints(vector<trial>& trial_vec) {
    trial_vec = trial_points;
}

void imgo_method::getPoints(vector<vector<double>> &points_vec) {
    points_vec.clear();
    vector<double> point(n);
    for (int i = 0; i < trial_points.size(); i++) {
        y(trial_points[i].x, point);
        points_vec.push_back(point);
    }
}

void imgo_method::y(double x, vector<double> &X) {
    mapd(x, den, X.data(), n, key);
    for (int i = 0; i < n; i++) {
        X[i] = X[i] * (B[i] - A[i]) + (A[i] + B[i]) / 2.0;
    }
}

double imgo_method::solve(int &n, Stop stop) {
    for (int i = 0; i < I.size(); i++) {
        I[i].clear();
    }
    trial_points.clear();
    M = 0;
    if (stop == NUMBER) Nmax = n;

    trial tr = newTrial(a);
    trial_points.push_back(tr);
    addInSort(I[tr.nu - 1], tr);
    tr = newTrial(b);
    trial_points.push_back(tr);
    addInSort(I[tr.nu - 1], tr);
    n = 2;

    double x_k_1;
    size_t t = 1;
    while(true) {
        x_k_1 = selectNewPoint(t);
        tr = newTrial(x_k_1);
        // Шаг 1
        addInSort(trial_points, tr);

        // Шаг 2
        addInSort(I[tr.nu - 1], tr);
        n++;
        if (stop == ACCURACY) {
            if (trial_points[t].x - trial_points[t - 1].x <= eps) {
                break;
            }
        } else {
            if (n >= Nmax) break;
        }

    }
    return searchMinX();
}

bool imgo_method::solve_test(double x_opt, int k) {
    for (int i = 0; i < I.size(); i++) {
        I[i].clear();
    }
    trial_points.clear();

    trial tr = newTrial(a);
    trial_points.push_back(tr);
    addInSort(I[tr.nu - 1], tr);
    tr = newTrial(b);
    trial_points.push_back(tr);
    addInSort(I[tr.nu - 1], tr);
    int n = 2;

    double x_k_1;
    size_t t = 1;
    while (true) {
        x_k_1 = selectNewPoint(t);
        tr = newTrial(x_k_1);
        // Шаг 1
        addInSort(trial_points, tr);

        // Шаг 2
        addInSort(I[tr.nu - 1], tr);
        n++;
        if (abs(x_k_1 - x_opt) <= eps) {
            return true;
        }
        if (n >= k) {
            return false;
        }
    }
}

void imgo_method::solve(int &count, vector<double> &X, Stop stop) {
    for (int i = 0; i < I.size(); i++) {
        I[i].clear();
    }
    trial_points.clear();
    M = 0;
    if (stop == NUMBER) Nmax = count;

    trial tr{0.0, -1.0, 0};
    trial_points.push_back(tr);
    tr.x = 1.0;
    trial_points.push_back(tr);
    tr = newTrial_md(0.5);
    addInSort(trial_points, tr);
    addInSort(I[tr.nu - 1], tr);
    if (tr.nu > M) {
        M = tr.nu;
    }
    count = 1;

    double x_k_1;
    size_t t = 1;
    while(true) {
        x_k_1 = selectNewPoint(t);
        tr = newTrial_md(x_k_1);
        // Шаг 1
        addInSort(trial_points, tr);

        // Шаг 2
        addInSort(I[tr.nu - 1], tr);
        if (tr.nu > M) {
            M = tr.nu;
        }

        count++;
        if (stop == ACCURACY) {
            if (pow(abs(trial_points[t].x - trial_points[t - 1].x), 1.0 / n) <= eps) {
                break;
            }
        } else {
            if (count >= Nmax) {
                break;
            }
        }
    }
    y(searchMinX(), X);
}
