#include <iostream>
#include <cmath>
#include <Grishagin/grishagin_function.hpp>
#include <GKLS/GKLSProblem.hpp>
#include <imgo.h>
#include <task.h>

using namespace std;

double f_test(vector<double> x, int j) {
    switch (j) {
        case 1: return 0.01 * (pow((x[0] - 2.2), 2.0) + pow((x[1] - 1.2), 2.0) - 2.25);
        case 2: return 100.0 * (1.0 - pow((x[0] - 2.0), 2.0) / 1.44 - pow(0.5 * x[1], 2.0));
        case 3: return 10.0 * (x[1] - 1.5 - 1.5 * sin(6.283 * (x[0] - 1.75)));
        case 4: return -1.5 * x[0] * x[0] * exp(1.0 - x[0] * x[0] - 20.25 * pow((x[0] - x[1]), 2.0)) - 
                       pow(0.5 * (x[0] - 1.0) * (x[1] - 1.0), 4.0) * exp(2.0 - pow(0.5 * (x[0] - 1.0), 4.0) - 
                       pow(x[1] - 1.0, 4.0));
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f_SPPR1(vector<double> x, int j) {
    switch (j) {
        case 1: return 1.0 - 2.0 * x[0] - x[0] * x[1];
        case 2: return 2.0 * x[0] * x[0] + x[1];
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f_SPPR2(vector<double> x, int j) {
    switch (j) {
        case 1: return sin(x[0]) - x[1];
        case 2: return cos(x[0]) + x[1];
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f_SPPR3(vector<double> x, int j) {
    switch (j) {
        case 1: return sin(x[0]) - x[1];
        case 2: return x[1] * cos(x[0]);
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f_SPPR4(vector<double> x, int j) {
    switch (j) {
        case 1: return 2 * sqrt(x[0]) - x[1];
        case 2: return -2 * x[0] * x[0] + x[1] * x[1];
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f_SPPR5(vector<double> x, int j) {
    switch (j) {
        case 1: return 0.5 * exp(x[0]) - x[1];
        case 2: return (x[0] - 2.0) * (x[0] - 2.0) + (log(x[1]) - 1) * (log(x[1]) - 1);
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f_SPPR6(vector<double> x, int j) {
    switch (j) {
        case 1: return 2.0 - sqrt(x[0]) - x[1];
        case 2: return exp(sqrt(x[0] * x[0] + (x[1] + 1.0) * (x[1] + 1.0)));
        default: return numeric_limits<double>::quiet_NaN();
    }
}

TGrishaginProblem grishaginProblem;
double f_grishagin(vector<double> x, int j) {
    switch (j) {
        case 1: return grishaginProblem.ComputeFunction(x);
        default: return numeric_limits<double>::quiet_NaN();
    }
}

TGKLSProblem gklsProblem1;
double f_gkls_1(vector<double> x, int j) {
    switch (j) {
        case 1: return gklsProblem1.ComputeFunction(x);
        default: return numeric_limits<double>::quiet_NaN();
    }
}

TGKLSProblem gklsProblem18(18);
double f_gkls_18(vector<double> x, int j) {
    switch (j) {
        case 1: return gklsProblem18.ComputeFunction(x);
        default: return numeric_limits<double>::quiet_NaN();
    }
}

int main() {
    int n = 2, m = 0, den = 10, key = 1, Nmax = 5000;
    double eps = 0.001, r = 2.0, d = 0.0;
    Stop stop = ACCURACY;
    vector<double> X(2);
    int number_trials;

    vector<double> A, B, A1, B1, A2, B2;
    grishaginProblem.GetBounds(A, B);
    gklsProblem1.GetBounds(A1, B1);
    gklsProblem18.GetBounds(A2, B2);
    vector<Task_peano> task{ Task_peano(f_test, "Test-1", n, 3, vector<double>{0.0, -1.0}, vector<double>{4.0, 3.0},
                                        vector<double>{0.942, 0.944}, eps, Nmax, r, d, den, key, stop, 0),
                             Task_peano(f_SPPR1, "SPPR-1", n, 1, vector<double>{0.5, -2.0}, vector<double>{3.0, 2.0},
                                        vector<double>{0.629960524947, -0.412598948032}, eps, Nmax, r, d, den, key, stop, 0),
                             Task_peano(f_SPPR2, "SPPR-2", n, 1, vector<double>{0.0, 0.0}, vector<double>{M_PI, 2.0},
                                        vector<double>{M_PI, 0.0}, eps, Nmax, 3.0, d, den, key, stop, 0),
                             Task_peano(f_SPPR3, "SPPR-3", n, 1, vector<double>{0.0, 0.0}, vector<double>{M_PI, 2.0},
                                        vector<double>{M_PI, 2.0}, eps, Nmax, r, d, den, key, stop, 0),
                             Task_peano(f_SPPR4, "SPPR-4", n, 1, vector<double>{0.0, 0.0}, vector<double>{4.0, 3.0},
                                        vector<double>{2.25, 3.0}, eps, Nmax, 2.5, d, den, key, stop, 0),
                             Task_peano(f_SPPR5, "SPPR-5", n, 1, vector<double>{0.0, 0.0}, vector<double>{10.0, 10.0},
                                        vector<double>{(3.0 + log(2.0)) / 2.0, 0.5 * exp((3.0 + log(2.0)) / 2.0)}, 
                                        eps, Nmax, 3.3, d, 12, key, ACCURNUMBER, 1),
                             Task_peano(f_SPPR6, "SPPR-6", n, 1, vector<double>{0.0, 0.0}, vector<double>{4.0, 2.0},
                                        vector<double>{1.0, 1.0}, eps, Nmax, 3.3, d, 12, key, stop, 1),
                             Task_peano(f_grishagin, "Grishagin_func", n, 0, A, B,
                                        grishaginProblem.GetOptimumPoint(), eps, Nmax, r, d, den, key, stop, 0),
                             Task_peano(f_gkls_1, "GKLS_func_1", n, 0, A1, B1,
                                        gklsProblem1.GetOptimumPoint(), eps, Nmax, 3.0, d, 12, key, stop, 0),
                             Task_peano(f_gkls_18, "GKLS_func_18", n, 0, A2, B2,
                                        gklsProblem18.GetOptimumPoint(), eps, Nmax, 4.0, d, den, key, stop, 0) };

    imgo_method imgo(nullptr, 2, 0, A, B);

    for (int i = 0; i < task.size(); i++) {
        if (task[i].used) {
            imgo.setFunc(task[i].f);
            imgo.setN(task[i].n);
            imgo.setM(task[i].m);
            imgo.setAB(task[i].A, task[i].B);
            imgo.setEps(task[i].eps);
            imgo.setNmax(task[i].Nmax);
            imgo.setR(task[i].r);
            imgo.setD(task[i].d);
            imgo.setDen(task[i].den);
            imgo.setKey(task[i].key);

            number_trials = task[i].Nmax;
            imgo.solve(number_trials, X, task[i].stop);

            cout << "Function: " << task[i].name << endl;
            cout << "Dimension = " << task[i].n << endl;
            cout << "Number of restrictions = " << task[i].m << endl;
            cout << "Parameters for constructing the Peano curve:" << endl;
            cout << "n = " << task[i].n << " m = " << task[i].den << " key = " << task[i].key << std::endl;
            cout << "Trials result:" << std::endl;
            cout << "Number of trials = " << number_trials << std::endl;
            cout << "x*_min = " << task[i].X_opt[0] << " y*_min = " << task[i].X_opt[1] << std::endl;
            cout << "x_min = " << X[0] << " y_min = " << X[1] << std::endl;
            cout << "||X* - X|| = " << sqrt((task[i].X_opt[0] - X[0]) * (task[i].X_opt[0] - X[0]) + 
                                            (task[i].X_opt[1] - X[1]) * (task[i].X_opt[1] - X[1])) << std::endl;
            cout << "|f(X*) - f(X)| = " << abs(task[i].f(task[i].X_opt, task[i].m + 1) - task[i].f(X, task[i].m + 1)) << std::endl;
            cout << std::endl;
        }
    }
    return 0;
}
