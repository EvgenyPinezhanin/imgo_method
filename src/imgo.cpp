#include<imgo.h>
#include<iostream>

imgo_method::imgo_method(double (*_f)(double, int), int _m, double _a, double _b, double _eps, double _r) 
    : f(_f), m(_m), a(_a), b(_b), eps(_eps), r(_r) {
    I.resize(m + 1);
    mu.resize(m + 1);
}

void imgo_method::addInSort(vector<trial> &vec, trial tr) {
    vector<trial>::iterator iter = vec.begin();
    vector<trial>::iterator iterEnd = vec.end();
    while(iter->x < tr.x) {
        if (iter != iterEnd) iter++;
    }
    vec.insert(iter, tr);
}

double imgo_method::searchMinX() {
    double z = a_inf[0].z;
    double x = a_inf[0].x;
    for (int i = 1; i < a_inf.size(); i++) {
        if (a_inf[i].z < z && a_inf[i].v == m + 1) {
            z = a_inf[i].z;
            x = a_inf[i].x;
        }
    }
    return x;
}

trial imgo_method::trial_func(double x) {
    int v = 1;
    double z = 0.0;
    for (int j = 1; j <= m + 1; j++) {
        if (v == m + 1) {
            z = f(x, v);
        }
        if (f(x, j) <= 0) {
            v++;
        } else {
            z = f(x, v);
        }
    }
    trial tr{x, z, v};
    return tr;
}

double imgo_method::selectNewPoint(int &t) {
    // Шаг 3
    double mu_tmp;
    for (int i = 0; i < m + 1; i++) {
        mu[i] = 0.0;
    }
    for (int i = 0; i < m + 1; i++) {
        if (I[i].size() < 2) {
            mu[i] = 1.0;
        } else {
            for (int j = 1; i < I[i].size(); j++) {
                mu_tmp = abs(I[i][j].z - I[i][j - 1].z) / (I[i][j].x - I[i][j - 1].x);
                if (mu_tmp > mu[i]) {
                    mu[i] = mu_tmp;
                }
            }
        }
    }
    for (int i = 0; i < m + 1; i++) {
        if (mu[i] == 0.0) {
            mu[i] = 1.0;
        };
    }

    // Шаг 4
    vector<double> z_v;

    double m = (M == 0.0) ? 1.0 : r * M;

    double tmp1 = a_inf[1].x - a_inf[0].x;
    double tmp2 = a_inf[1].z - a_inf[0].z;
    double tmp3 = a_inf[1].z + a_inf[0].z;
    double R = m * tmp1 + (tmp2 * tmp2) / (m * tmp1) - 2 * tmp3;
    double Rtmp;
    t = 1;
    int n = a_inf.size();
    for (int i = 2; i < n; i++) {
        tmp1 = a_inf[i].x - a_inf[i - 1].x;
        tmp2 = a_inf[i].z - a_inf[i - 1].z;
        tmp3 = a_inf[i].z + a_inf[i - 1].z;
        Rtmp = m * tmp1 + (tmp2 * tmp2) / (m * tmp1) - 2 * tmp3;
        if (Rtmp > R) {
            R = Rtmp;
            t = i;
        }
    }
    x_k_1 = (a_inf[t].x + a_inf[t - 1].x) / 2 - (a_inf[t].z - a_inf[t - 1].z) / (2 * m);
    if (x_k_1 < a_inf[t - 1].x || x_k_1 > a_inf[t].x) {
        cout << x_k_1;
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

void imgo_method::setM(int _m) {
    m = _m;
    if (I.size() < m + 1) {
        I.resize(m + 1);
    }
}

double imgo_method::solve(int &n) {
    for (int i = 0; I.size(); i++) {
        I[i].clear();
    }
    a_inf.clear();
    a_inf.push_back(trial_func(a));
    I[a_inf[0].v - 1].push_back(a_inf[0]);
    a_inf.push_back(trial_func(b));
    I[a_inf[1].v - 1].push_back(a_inf[1]);
    n = 2;

    double x_k_1;
    int t = 1;
    trial tr = {0.0, 0.0, 0};
    while(true) {
        tr = trial_func(selectNewPoint(t));
        // Шаг 1
        addInSort(a_inf, tr);
        // Шаг 2
        addInSort(I[tr.v], tr);
        n++;
        if (a_inf[t].x - a_inf[t - 1].x <= eps) {
            break;
        }
    }
    return searchMinX();
}