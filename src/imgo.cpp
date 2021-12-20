#include<imgo.h>
#include<iostream>

imgo_method::imgo_method(double (*_f)(double, int), size_t _m, double _a, double _b, double _eps, double _r, double _d) 
    : f(_f), m(_m), a(_a), b(_b), eps(_eps), r(_r), d(_d) {
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
                mu_tmp = abs(I[nu][i].z - I[nu][j].z) / (I[nu][i].x - I[nu][j].x);
                if (mu_tmp > mu[nu]) {
                    mu[nu] = mu_tmp;
                }
            }
        }
        if (abs(mu[nu]) < 6e-13) {
            mu[nu] = 1.0;
        };
    }

    // Шаг 4
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

    // Шаг 5, 6
    t = 1;
    double d_x = trial_points[1].x - trial_points[0].x; 
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
    size_t size_tr_pt = trial_points.size();
    for (size_t i = 2; i < size_tr_pt; i++) {
        d_x = trial_points[i].x - trial_points[i - 1].x;
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
    }

    // Шаг 7
    double x_k_1;
    if (trial_points[t].nu != trial_points[t - 1].nu) {
        x_k_1 = (trial_points[t].x + trial_points[t - 1].x) / 2.0;
    } else {
        x_k_1 = (trial_points[t].x + trial_points[t - 1].x) / 2.0 
                - (trial_points[t].z - trial_points[t - 1].z) / (2.0 * r * mu[trial_points[t].nu - 1]);
    }
    return x_k_1;
}

void imgo_method::setFunc(double (*_f)(double, int)) {
    f = _f;
}

void imgo_method::setA(double _a) {
    a = _a;
}

void imgo_method::setB(double _b) {
    b = _b;
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

void imgo_method::getTrialPoints(vector<trial>& trial_vec) {
    trial_vec = trial_points;
}

double imgo_method::solve(int &n) {
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
        if (trial_points[t].x - trial_points[t - 1].x <= eps) {
            break;
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