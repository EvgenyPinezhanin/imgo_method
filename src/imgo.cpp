#include<imgo.h>
#include<iostream>

imgo_method::imgo_method(double (*_f)(double, int), int _m, double _a, double _b, double _eps, double _r) 
    : f(_f), m(_m), a(_a), b(_b), eps(_eps), r(_r) {
    I.resize(m + 1);
    mu.resize(m + 1);
    z_v.resize(m + 1);
}

void imgo_method::addInSort(vector<trial> &vec, trial tr) {
    if (vec.empty()) {
        vec.push_back(tr);
        return;
    }
    vector<trial>::iterator iter = vec.begin();
    vector<trial>::iterator iterEnd = vec.end();
    while(iter->x < tr.x) {
        iter++;
        if (iter == iterEnd) break;
    }
    vec.insert(iter, tr);
}

double imgo_method::searchMinX() {
    double z;
    double x;
    int k;
    for (int i = 1; i < a_inf.size(); i++) {
        if (a_inf[i].v == m + 1) {
            z = a_inf[i].z;
            x = a_inf[i].x;
            k = i + 1;
            break;
        }
    }
    for (int i = k; i < a_inf.size(); i++) {
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
            break;
        }
        if (f(x, j) <= 0) {
            v++;
        } else {
            z = f(x, v);
            break;
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
    for (int v = 0; v < m + 1; v++) {
        for (int i = 1; i < I[v].size(); i++) {
            for (int j = 0; j < i; j++){
                mu_tmp = abs(I[v][i].z - I[v][j].z) / (I[v][i].x - I[v][j].x);
                if (mu_tmp > mu[v]) {
                    mu[v] = mu_tmp;
                }
            }
        }
    }
    for (int v = 0; v < m + 1; v++) {
        if (mu[v] == 0.0) {
            mu[v] = 1.0;
        };
    }

    // Шаг 4
    for (int v = 0; v < m + 1; v++) {
        if (I[v].size() != 0) {
            z_v[v] = I[v][0].z;
            for (int i = 1; i < I[v].size(); i++) {
                if (I[v][i].z < z_v[v]) {
                    z_v[v] = I[v][i].z;
                }
            }
        }
        if (z_v[v] <= 0 && v != m) {
            z_v[v] = -eps;
        }
    }

    // Шаг 5, 6
    t = 1;
    double d_x = a_inf[1].x - a_inf[0].x; 
    double R, Rtmp;
    double mu_v, z_v_v;
    if (a_inf[1].v == a_inf[0].v) {
        mu_v = mu[a_inf[1].v - 1];
        z_v_v = z_v[a_inf[1].v - 1];
        R = d_x + pow(a_inf[1].z - a_inf[0].z, 2) / (r * r * mu_v * mu_v * d_x) -
        2.0 * (a_inf[1].z + a_inf[0].z - 2.0 * z_v_v) / (r * mu_v);
    } else if (a_inf[0].v < a_inf[1].v) {
        mu_v = mu[a_inf[1].v - 1];
        z_v_v = z_v[a_inf[1].v - 1];
        R = 2 * d_x  - 4.0 * (a_inf[1].z - z_v_v) / (r * mu_v);
    } else if (a_inf[0].v > a_inf[1].v) {
        mu_v = mu[a_inf[0].v - 1];
        z_v_v = z_v[a_inf[0].v - 1];
        R = 2 * d_x  - 4.0 * (a_inf[0].z - z_v_v) / (r * mu_v);
    }
    int n = a_inf.size();
    for (int i = 2; i < n; i++) {
        d_x = a_inf[i].x - a_inf[i - 1].x; 
        if (a_inf[i].v == a_inf[i - 1].v) {
            mu_v = mu[a_inf[i].v - 1];
            z_v_v = z_v[a_inf[i].v - 1];
            Rtmp = d_x + pow(a_inf[i].z - a_inf[i - 1].z, 2) / (r * r * mu_v * mu_v * d_x) -
            2.0 * (a_inf[i].z + a_inf[i - 1].z - 2.0 * z_v_v) / (r * mu_v);
        } else if (a_inf[i - 1].v < a_inf[i].v) {
            mu_v = mu[a_inf[i].v - 1];
            z_v_v = z_v[a_inf[i].v - 1];
            Rtmp = 2 * d_x  - 4.0 * (a_inf[i].z - z_v_v) / (r * mu_v);
        } else if (a_inf[i - 1].v > a_inf[i].v) {
            mu_v = mu[a_inf[i - 1].v - 1];
            z_v_v = z_v[a_inf[i - 1].v - 1];
            Rtmp = 2 * d_x  - 4.0 * (a_inf[i - 1].z - z_v_v) / (r * mu_v);
        }
        if (Rtmp > R) {
            R = Rtmp;
            t = i;
        }
    }

    // Шаг 7
    double x_k_1;
    if (a_inf[t].v != a_inf[t - 1].v) {
        x_k_1 = (a_inf[t].x + a_inf[t - 1].x) / 2.0;
    } else {
        x_k_1 = (a_inf[t].x + a_inf[t - 1].x) / 2.0 - (a_inf[t].z - a_inf[t - 1].z) / (2 * r * mu[a_inf[t].v - 1]);
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
    I.resize(m + 1);
    mu.resize(m + 1);
    z_v.resize(m + 1);
}

double imgo_method::solve(int &n) {
    for (int i = 0; i < I.size(); i++) {
        I[i].clear();
    }
    a_inf.clear();

    trial tr = trial_func(a);
    a_inf.push_back(tr);
    I[tr.v - 1].push_back(tr);
    tr = trial_func(b);
    a_inf.push_back(tr);
    I[tr.v - 1].push_back(tr);
    n = 2;

    double x_k_1;
    int t = 1;
    while(true) {
        x_k_1 = selectNewPoint(t);
        tr = trial_func(x_k_1);
        // Шаг 1
        addInSort(a_inf, tr);
        // Шаг 2
        for (int i = 0; i < tr.v; i++) {
            addInSort(I[i], tr);
        }
        n++;
        if (a_inf[t].x - a_inf[t - 1].x <= eps) {
            break;
        }
    }
    return searchMinX();
}