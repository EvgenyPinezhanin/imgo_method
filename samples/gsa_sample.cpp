#if defined( _MSC_VER )
    #define _CRT_SECURE_NO_WARNINGS
#endif

#include <iostream>
#include <fstream>
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

const int sample_func_number = 3; // 0 - f1, 1 - f2, 2 - f3, 3 - f4

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
    ofstream ofstr("output_data/gsa_sample.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";

    double x, eps = 0.0001, r = 2.0;
    int count_trials, Nmax = 1000;
    Stop stop = Stop::ACCURACY;

    vector<task_gsa> task_array = { task_gsa(f1, "f1(x) = -4.0 * x + 1.0", 3.0, 4.0, 4.0, 4.0, eps, Nmax, r, stop),
                                    task_gsa(f2, "f2(x) = 5.0 * x * x + 3.0 * x - 1.0", -2.0, 2.0, -0.3, 23.0, eps, Nmax, r, stop),
                                    task_gsa(f3, "f3(x) = x * sin(x)", 0.0, 20.0, 17.336, 18.955, eps ,Nmax, 2.1, stop),
                                    task_gsa(f4, "f4(x) = x * sin(1 / x)", -0.4, -0.05, -0.2225, 6 * M_PI, eps, Nmax, r, stop) };

    gsa_method gsa(nullptr);

    vector<trial> trials;
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
            cout << "Lipschitz constant = " << task_array[i].L[0] << endl;
            cout << "X* = " << setprecision(8) << task_array[i].X_opt[0] << endl;
            cout << "f(X*) = " << setprecision(8) << task_array[i].f(task_array[i].X_opt[0]) << endl;
            cout << "Parameters for method:" << endl;
            cout << "eps = " << eps << " r = " << r << endl;
            cout << "Trials result:" << endl;
            cout << "Number of trials = " << count_trials << endl;
            cout << "Estimation of the Lipschitz constant = " << gsa.getLambda() << endl;
            cout << "X = " << setprecision(8) << x << endl;
            cout << "f(X) = " << setprecision(8) << task_array[i].f(x) << endl;
            cout << "|X* - X| = " << setprecision(8) << abs(task_array[i].X_opt[0] - x) << endl;
            cout << "|f(X*) - f(X)| = " << setprecision(8) << abs(task_array[i].f(task_array[i].X_opt[0]) - task_array[i].f(x)) << endl;
            cout << endl;

            // Saving points for plotting
            ofstr << task_array[i].X_opt[0] << " " << task_array[i].f(task_array[i].X_opt[0]) << endl;
            ofstr << endl << endl;
            ofstr << x << " " << task_array[i].f(x) << endl;
            ofstr << endl << endl;
            gsa.getTrialPoints(trials);
            for (int j = 0; j < trials.size(); j++) {
                ofstr << trials[j].x << " " << trials[j].z << endl;
            }
            ofstr << endl << endl;
        }
    }
    ofstr.close();

    // Plotting the function(works with gnuplot)
    int error;
#if defined(__linux__)
    setenv("QT_QPA_PLATFORM", "xcb", false);
    error = system("chmod +x scripts/gsa_sample.gp");
    if (error != 0) {
        cerr << "Error chmod" << endl;
    }
#endif

    char str[100];
    sprintf(str, "gnuplot -c scripts/gsa_sample.gp %d", sample_func_number);
    error = system(str);
    if (error != 0) {
        cerr << "Error gnuplot" << endl;
    }

#if defined( _MSC_VER )
    cin.get();
#endif

    return 0;
}
