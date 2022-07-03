#include<gsa.h>

gsa_method::gsa_method(double (*_f)(double), double _a, double _b, double _eps, double _r) 
    : f(_f), a(_a), b(_b), eps(_eps), r(_r) { }

void gsa_method::addInSort(double x) {
    vector<trial>::iterator iter = trial_points.begin();
    vector<trial>::iterator iterEnd = trial_points.end();
    while (true) {
        if (iter == iterEnd || iter->x > x) break;
        iter++;
    }
    trial_points.insert(iter, trial(x, f(x)));
}

double gsa_method::searchMinX() {
    double z = trial_points[0].z;
    double x = trial_points[0].x;
    for (int i = 1; i < trial_points.size(); i++) {
        if (trial_points[i].z < z) {
            z = trial_points[i].z;
            x = trial_points[i].x;
        }
    }
    return x;
}

double gsa_method::selectNewPoint(size_t &t) {
    static double x_k_1 = b;
    static double M = -1.0;
    double tmpM = 0.0;

    // Шаг 2
    if (trial_points[t].x == b) {
        M = max({M, abs((f(x_k_1) - trial_points[t - 1].z) / (x_k_1 - trial_points[t - 1].x))});
    } else {
        M = max({M, abs((f(x_k_1) - trial_points[t - 1].z) / (x_k_1 - trial_points[t - 1].x)), abs((trial_points[t + 1].z - f(x_k_1)) / (trial_points[t + 1].x - x_k_1))});
    }

    // Шаг 3
    double m = (M == 0.0) ? 1.0 : r * M;

    // Шаг 4, 5 
    double d_x = trial_points[1].x - trial_points[0].x;
    double R = m * d_x + pow(trial_points[1].z - trial_points[0].z, 2.0) / (m * d_x) - 2 * (trial_points[1].z + trial_points[0].z);
    double Rtmp;
    t = 1;
    size_t n = trial_points.size();
    for (size_t i = 2; i < n; i++) {
        d_x = trial_points[i].x - trial_points[i - 1].x;
        Rtmp = m * d_x + pow(trial_points[i].z - trial_points[i - 1].z, 2.0) / (m * d_x) - 2 * (trial_points[i].z + trial_points[i - 1].z);
        if (Rtmp > R) {
            R = Rtmp;
            t = i;
        }
    }

    // Шаг 6
    x_k_1 = (trial_points[t].x + trial_points[t - 1].x) / 2 - (trial_points[t].z - trial_points[t - 1].z) / (2 * m);
    return x_k_1;
}

void gsa_method::setFunc(double (*_f)(double)) {
    f = _f;
}

void gsa_method::setA(double _a) {
    a = _a;
}

void gsa_method::setB(double _b) {
    b = _b;
}

void gsa_method::setEps(double _eps) {
    eps = _eps;
}

double gsa_method::solve(int &n) {
    trial_points.clear();

    trial_points.push_back(trial(a, f(a)));
    trial_points.push_back(trial(b, f(b)));
    n = 2;

    double x_k_1;
    size_t t = 1;
    while(true) {
        x_k_1 = selectNewPoint(t);

        // Шаг 1
        addInSort(x_k_1);
        n++;
        if (trial_points[t].x - trial_points[t - 1].x <= eps) {
            break;
        }
    }
    return searchMinX();
}