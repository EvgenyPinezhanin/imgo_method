#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#if defined( _MSC_VER )
    #define _USE_MATH_DEFINES
    #include <math.h>
#else
    #include <cmath>
#endif

#include <imgo.h>
#include <task.h>

using namespace std;

double f1(double x, int j) {
    switch(j) {
        case 1: return sin(x);
        case 2: return -2.0 * x + 3.0;
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f2(double x, int j) {
    switch (j) {
        case 1: return x * x - 0.05;
        case 2: return -x + 0.1;
        case 3: return 5.0 * x * x + 3.0 * x - 1.0;
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f3(double x, int j) {
    switch (j) {
        case 1: return sin(x);
        default: return numeric_limits<double>::quiet_NaN();
    }
}

int main() {
    double x, eps = 0.001, r = 3.0, d = 0.0;
    int count_trials, Nmax = 1000;
    Stop stop = Stop::ACCURACY;

    vector<task_imgo> task_array = { task_imgo(f1, "f1(x)", 1, 2.0, 8.0, 2.0 * M_PI, eps, Nmax, r, d, stop),
                                     task_imgo(f2, "f2(x)", 2, -2.0, 2.0, 0.1, eps, Nmax, r, d, stop),
                                     task_imgo(f3, "f3(x)", 0, -4.0, 4.0, -M_PI / 2.0, eps ,Nmax, r, d, stop) };

    imgo_method imgo(nullptr);

    vector<double> mu;
    for (int i = 0; i < task_array.size(); i++) {
        imgo.setF(task_array[i].f);
        imgo.setM(task_array[i].m);
        imgo.setAB(task_array[i].A[0], task_array[i].B[0]);
        imgo.setEps(task_array[i].eps);
        imgo.setNmax(task_array[i].Nmax);
        imgo.setR(task_array[i].r);
        imgo.setD(task_array[i].d);

        imgo.solve(count_trials, x, task_array[i].stop);
        imgo.getMu(mu);

        cout << "Function: " << task_array[i].name << endl;
        cout << "Number of constrained = " << task_array[i].m << endl;
        cout << "[a; b] = [" << task_array[i].A[0] << "; " << task_array[i].B[0] << "]"<< endl;
        cout << "X* = " << setprecision(8) << task_array[i].X_opt[0] << endl;
        cout << "f(X*) = " << setprecision(8) << task_array[i].f(task_array[i].X_opt[0], task_array[i].m + 1) << endl;
        cout << "Parameters for method:" << endl;
        cout << "eps = " << eps << " r = " << r << " d = " << d << endl;
        cout << "Trials result:" << endl;
        cout << "Number of trials = " << count_trials << endl;
        cout << "Estimation of the Lipschitz constant:" << endl;
        cout << "L(" << task_array[i].name << ") = " << mu[task_array[i].m] << endl;
        for (int j = 0; j < task_array[i].m; j++) {
            cout << "L(g" << j + 1 << ") = " << mu[j] << endl;
        }
        cout << "X = " << setprecision(8) << x << endl;
        cout << "f(X) = " << setprecision(8) << task_array[i].f(x, task_array[i].m + 1) << endl;
        cout << "|X* - X| = " << setprecision(8) << abs(task_array[i].X_opt[0] - x) << endl;
        cout << "|f(X*) - f(X)| = " << setprecision(8) << abs(task_array[i].f(task_array[i].X_opt[0], task_array[i].m + 1) - 
                                                              task_array[i].f(x, task_array[i].m + 1)) << endl;
        cout << endl;
    }

#if defined( _MSC_VER )
    cin.get();
#endif

    return 0;
}
