#include <iostream>
#include <iomanip>
#include <vector>
#if defined( _MSC_VER )
    #define _USE_MATH_DEFINES
    #include <math.h>
#else
    #include <cmath>
#endif

#include <gsa.h>
#include <task.h>

using namespace std;

double f1(double x) {
    return -4.0 * x + 1.0;
}

double f2(double x) {
    return 5.0 * x * x + 3.0 * x - 1.0;
}

double f3(double x) {
    return x * sin(x);
}

double f4(double x) {
    if (x != 0) {
        return x * sin(1 / x);
    } else {
        return 0.0;
    }
}

int main() {
    double x, eps = 0.0001, r = 2.1;
    int count_trials, Nmax = 1000;
    Stop stop = Stop::ACCURACY;

    vector<task_gsa> task_array = { task_gsa(f1, "f1(x) = -4.0 * x + 1.0", 3.0, 4.0, 4.0, eps, Nmax, r, stop),
                                    task_gsa(f2, "f2(x) = 5.0 * x * x + 3.0 * x - 1.0", -2.0, 2.0, -0.3, eps, Nmax, r, stop),
                                    task_gsa(f3, "f3(x) = x * sin(x)", 0.0, 20.0, 17.336, eps ,Nmax, r, stop),
                                    task_gsa(f4, "f4(x) = x * sin(1 / x)", -0.4, 0.4, -0.2225, eps, Nmax, r, stop) };

    gsa_method gsa(nullptr);

    for (int i = 0; i < task_array.size(); i++) {
        if (task_array[i].used) {
            gsa.setF(task_array[i].f);
            gsa.setAB(task_array[i].A[0], task_array[i].B[0]);
            gsa.setEps(task_array[i].eps);
            gsa.setNmax(task_array[i].Nmax);
            gsa.setR(task_array[i].r);

            gsa.solve(count_trials, x, task_array[i].stop);

            cout << "Function: " << task_array[i].name << endl;
            cout << "[a; b] = [" << task_array[i].A[0] << "; " << task_array[i].B[0] << "]"<< endl;
            cout << "X* = " << setprecision(8) << task_array[i].X_opt[0] << endl;
            cout << "f(X*) = " << setprecision(8) << task_array[i].f(task_array[i].X_opt[0]) << endl;
            cout << "Parameters for method:" << endl;
            cout << "eps = " << eps << " r = " << r << endl;
            cout << "Trials result:" << endl;
            cout << "Number of trials = " << count_trials << endl;
            cout << "Estimation of the Lipschitz constant = " << gsa.getMu() << endl;
            cout << "X = " << setprecision(8) << x << endl;
            cout << "f(X) = " << setprecision(8) << task_array[i].f(x) << endl;
            cout << "|X* - X| = " << setprecision(8) << abs(task_array[i].X_opt[0] - x) << endl;
            cout << "|f(X*) - f(X)| = " << setprecision(8) << abs(task_array[i].f(task_array[i].X_opt[0]) - task_array[i].f(x)) << endl;
            cout << endl;
        }
    }

#if defined( _MSC_VER )
    cin.get();
#endif

    return 0;
}
