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
        //cout << "Mu[nu] = " << mu[nu] << endl;
        if (abs(mu[nu]) < 6e-13) {
            mu[nu] = 1.0;
        };
        //cout << "Mu[nu] = " << mu[nu] << endl;
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
            //cout << "z_star[nu] = " << z_star[nu] << " nu = " << nu << endl;
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

//test
void imgo_method::getTrialPointsTest(vector<trial>& trial_vec) {
    trial_vec = trial_points_test;
}
//test

double imgo_method::solve(int &n) {
    for (int i = 0; i < I.size(); i++) {
        I[i].clear();
    }
    trial_points.clear();
    trial_points_test.clear();

    trial tr = newTrial(a);
    trial_points.push_back(tr);
    //test
    trial_points_test.push_back(tr);
    //test
    addInSort(I[tr.nu - 1], tr);
    tr = newTrial(b);
    trial_points.push_back(tr);
    //test
    trial_points_test.push_back(tr);
    //test
    addInSort(I[tr.nu - 1], tr);
    n = 2;

    double x_k_1;
    size_t t = 1;
    while(true) {
        x_k_1 = selectNewPoint(t);
        tr = newTrial(x_k_1);
        // Шаг 1
        addInSort(trial_points, tr);

        //test
        trial_points_test.push_back(tr);
        //test

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
    trial_points_test.clear();

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




/*
imgo_method_adaptive::imgo_method_adaptive(double (*_f)(double, int), int _m, double _a, double _b, double _eps, double _r, double _d) 
    : imgo_method(_f, _m, _a, _b, _eps, _r, _d) { 
    H_default.resize(m + 1);
    for (int i = 1; i <= m + 1; i++) {
        H_default[i - 1] = i;
    }
}

void imgo_method_adaptive::changeH(vector<int> &H, int p) {
    for (int i = p - 2; i >= 0; i--) {
        H[i + 1] = H[i];
    }
    H[0] = p;
}

void imgo_method_adaptive::changeH(vector<int> &H, int p, int q) {
    int pt = p, qt = q;
    if (pt > qt) swap(pt, qt);
    for (int i =  qt - 2; i > pt - 1; i--) {
        H[i + 1] = H[i];
    }
    for (int i =  pt - 2; i >= 0; i--) {
        H[i + 2] = H[i];
    }
    H[1] = p;
    H[0] = q;
}

// от 1 до m+1
int imgo_method_adaptive::getNumberOfIndex(double v, trial tr) {
    for(int i = 0; i < m + 1; i++) {
        if (v == H[tr.H_number][i]) {
            return i + 1;
        }
    }
    return -1;
}

void imgo_method_adaptive::trial_func(trial &tr) {
    int v = 1;
    double z = 0.0;
    for (int j = 1; j <= m + 1; j++) {
        if (j == m + 1) {
            z = f(tr.x, H[tr.H_number][j - 1]);
            v = H[tr.H_number][j - 1];
            break;
        }
        if (f(tr.x, j) <= 0) {

        } else {
            z = f(tr.x, H[tr.H_number][j - 1]);
            v = H[tr.H_number][j - 1];
            break;
        }
    }
    tr.z = z;
    tr.v = v;
}

trial imgo_method_adaptive::selectNewPoint(int &t) {
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
            z_v[v] = -mu[v] * d;
        }
    }

    // Шаг 5, 6
    t = 1;
    int p, q, u, s;
    double d_x = a_inf[1].x - a_inf[0].x; 
    double R, Rtmp;
    double mu_v, z_v_v;
    if (a_inf[0].v == a_inf[1].v) {
        mu_v = mu[a_inf[1].v - 1];
        z_v_v = z_v[a_inf[1].v - 1];
        R = d_x + pow(a_inf[1].z - a_inf[0].z, 2) / (r * r * mu_v * mu_v * d_x) -
        2.0 * (a_inf[1].z + a_inf[0].z - 2.0 * z_v_v) / (r * mu_v);
    } else {
        p = getNumberOfIndex(a_inf[0].v, a_inf[0]);
        q = getNumberOfIndex(a_inf[0].v, a_inf[1]);
        u = getNumberOfIndex(a_inf[1].v, a_inf[0]);
        s = getNumberOfIndex(a_inf[1].v, a_inf[1]);
        if (p < u && q < s) {
            mu_v = mu[a_inf[1].v - 1];
            z_v_v = z_v[a_inf[1].v - 1];
            R = 2 * d_x  - 4.0 * (a_inf[1].z - z_v_v) / (r * mu_v);
        } else if (p > u && q > s) {
            mu_v = mu[a_inf[0].v - 1];
            z_v_v = z_v[a_inf[0].v - 1];
            R = 2 * d_x  - 4.0 * (a_inf[0].z - z_v_v) / (r * mu_v);
        } else {
            if (a_inf[0].v < a_inf[1].v) {
                mu_v = mu[a_inf[1].v - 1];
                z_v_v = z_v[a_inf[1].v - 1];
                R = 2 * d_x  - 4.0 * (a_inf[1].z - z_v_v) / (r * mu_v);
            } else if (a_inf[0].v > a_inf[1].v) {
                mu_v = mu[a_inf[0].v - 1];
                z_v_v = z_v[a_inf[0].v - 1];
                R = 2 * d_x  - 4.0 * (a_inf[0].z - z_v_v) / (r * mu_v);
            }
        }
    }

    int n = a_inf.size();
    for (int i = 2; i < n; i++) {
        d_x = a_inf[i].x - a_inf[i - 1].x; 
        if (a_inf[i].v == a_inf[i - 1].v) {
            mu_v = mu[a_inf[i].v - 1];
            z_v_v = z_v[a_inf[i].v - 1];
            Rtmp = d_x + pow(a_inf[i].z - a_inf[i - 1].z, 2) / (r * r * mu_v * mu_v * d_x) -
            2.0 * (a_inf[i].z + a_inf[i - 1].z - 2.0 * z_v_v) / (r * mu_v);
        } else {
            p = getNumberOfIndex(a_inf[i - 1].v, a_inf[i - 1]);
            q = getNumberOfIndex(a_inf[i - 1].v, a_inf[i]);
            u = getNumberOfIndex(a_inf[i].v, a_inf[i - 1]);
            s = getNumberOfIndex(a_inf[i].v, a_inf[i]);
            if (p < u && q < s) {
                mu_v = mu[a_inf[i].v - 1];
                z_v_v = z_v[a_inf[i].v - 1];
                Rtmp = 2 * d_x  - 4.0 * (a_inf[i].z - z_v_v) / (r * mu_v);
            } else if (p > u && q > s) {
                mu_v = mu[a_inf[i - 1].v - 1];
                z_v_v = z_v[a_inf[i - 1].v - 1];
                Rtmp = 2 * d_x  - 4.0 * (a_inf[i - 1].z - z_v_v) / (r * mu_v);
            } else {
                if (a_inf[i - 1].v < a_inf[i].v) {
                    mu_v = mu[a_inf[i].v - 1];
                    z_v_v = z_v[a_inf[i].v - 1];
                    Rtmp = 2 * d_x  - 4.0 * (a_inf[i].z - z_v_v) / (r * mu_v);
                } else if (a_inf[i - 1].v > a_inf[i].v) {
                    mu_v = mu[a_inf[i - 1].v - 1];
                    z_v_v = z_v[a_inf[i - 1].v - 1];
                    Rtmp = 2 * d_x  - 4.0 * (a_inf[i - 1].z - z_v_v) / (r * mu_v);
                }
            }
        }
        if (Rtmp > R) {
            R = Rtmp;
            t = i;
        }
    }

    // Шаг 7
    trial tr;
    double x_k_1;
    if (a_inf[t].v != a_inf[t - 1].v) {
        x_k_1 = (a_inf[t].x + a_inf[t - 1].x) / 2.0;
    } else {
        x_k_1 = (a_inf[t].x + a_inf[t - 1].x) / 2.0 - (a_inf[t].z - a_inf[t - 1].z) / (2 * r * mu[a_inf[t].v - 1]);
    }
    tr.x = x_k_1;
    p = a_inf[t - 1].v;
    q = a_inf[t].v;
    vector<int> H_new(H_default);
    if (p == q && p == m + 1) {

    } else if ((p == q && p < m + 1) || (p < q && p == m + 1)) {
        changeH(H_new, p);
    } else if (q < p && p == m + 1) {
        changeH(H_new, q);
    } else if (p < q && q <= m) {
        changeH(H_new, q, p);
    } else if (q < p && p <= m) {
        changeH(H_new, p, q);
    }
    H.push_back(H_new);
    // Debug
    //if (H_new[0] == 1 && H_new[1] == 2) {
    //
    //} else {
    //    cout << "Error" << endl;
    //}
    //
    tr.H_number = H.size() - 1;
    return tr;
}

void imgo_method_adaptive::setM(int _m) {
    imgo_method::setM(_m);
    H_default.resize(m + 1);
    for (int i = 1; i <= m + 1; i++) {
        H_default[i - 1] = i;
    }
}

double imgo_method_adaptive::solve(int &n) {
    for (int i = 0; i < I.size(); i++) {
        I[i].clear();
    }
    a_inf.clear();
    H.clear();

    trial tr(a);
    H.push_back(H_default);
    tr.H_number = H.size() - 1;
    trial_func(tr);
    a_inf.push_back(tr);
    I[(size_t)tr.v - 1].push_back(tr);
    H.push_back(H_default);
    tr.H_number = H.size() - 1;
    tr.x = b;
    trial_func(tr);
    a_inf.push_back(tr);
    I[(size_t)tr.v - 1].push_back(tr);
    n = 2;

    int t = 1;
    while(true) {
        tr = selectNewPoint(t);
        trial_func(tr);
        // Шаг 1
        addInSort(a_inf, tr);
        // Шаг 2
        for (int i = 0; i < tr.v; i++) {
            addInSort(I[i], tr);
        }
        n++;
        if (a_inf[t].x - a_inf[(size_t)t - 1].x <= eps) {
            break;
        }
    }
    return searchMinX();
}*/