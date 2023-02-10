#include <iostream>
#include <fstream>
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

const int sample_func_number = 2; // 0 - f1, 1 - f2, 2 - f3

double f1(double x, int j) {
    switch (j) {
        case 1: return sin(x);
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f2(double x, int j) {
    switch(j) {
        case 1: return sin(x);
        case 2: return -2.0 * x + 3.0;
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f3(double x, int j) {
    switch (j) {
        case 1: return x * x - 0.05;
        case 2: return -x + 0.1;
        case 3: return 5.0 * x * x + 3.0 * x - 1.0;
        default: return numeric_limits<double>::quiet_NaN();
    }
}

int main() {
    ofstream ofstr("output_data/imgo_sample.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstr_opt("output_data/imgo_sample_opt.txt");
    if (!ofstr_opt.is_open()) cerr << "File opening error\n";

    double x, eps = 0.001, r = 3.0, d = 0.0;
    int count_trials, Nmax = 1000;
    Stop stop = Stop::ACCURACY;

    vector<task_imgo> task_array = { task_imgo(f1, "f1(x)", 0, -4.0, 4.0, -M_PI / 2.0, vector<double>{1.0}, eps ,Nmax, r, d, stop),
                                     task_imgo(f2, "f2(x)", 1, 2.0, 8.0, 2.0 * M_PI, vector<double>{1.0 ,2.0}, eps, Nmax, r, d, stop),
                                     task_imgo(f3, "f3(x)", 2, -2.0, 2.0, 0.1, vector<double>{4.0, 1.0, 23.0}, eps, Nmax, r, d, stop) };

    imgo_method imgo(nullptr);

    vector<double> mu;
    vector<trial_constr> trials;
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
        cout << "Lipschitz constant:" << endl;
        cout << "L*(" << task_array[i].name << ") = " << task_array[i].L[task_array[i].m] << endl;
        for (int j = 0; j < task_array[i].m; j++) {
            cout << "L*(g" << j + 1 << ") = " << task_array[i].L[j] << endl;
        }
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

        // Saving points for plotting
        ofstr << task_array[i].X_opt[0] << " " << task_array[i].f(task_array[i].X_opt[0], task_array[i].m + 1) << endl;
        ofstr << endl << endl;
        ofstr << x << " " << task_array[i].f(x, task_array[i].m + 1) << endl;
        ofstr << endl << endl;
        imgo.getTrialPoints(trials);
        for (int j = 0; j < trials.size(); j++) {
            ofstr << trials[j].x << " " << trials[j].z << endl;
        }
        ofstr << endl << endl;
    }
    ofstr.close();

    ofstr_opt << "array functions[" << task_array.size() << "]" << endl;
    for (int i = 0; i < task_array.size(); i++) {
        ofstr_opt << "functions[" << i + 1 << "] = \"f" << i + 1 << "(x) title \\\"f(x)\\\"";
        for (int j = 0; j < task_array[i].m; j++) {
            ofstr_opt << ", g" << i + 1 << "_" << j + 1 << "(x) title \\\"g" << j + 1 <<"(x)\\\"";
        }
        ofstr_opt << "\""<< endl;
    }
    ofstr_opt.close();

    // Plotting the function(works with gnuplot)
    int error;
#if defined(__linux__)
    setenv("QT_QPA_PLATFORM", "xcb", false);
    error = system("chmod +x scripts/imgo_sample.gp");
    if (error != 0) {
        cerr << "Error chmod" << endl;
    }
#endif

    char str[100];
    sprintf(str, "gnuplot -c scripts/imgo_sample.gp %d", sample_func_number);
    error = system(str);
    if (error != 0) {
        cerr << "Error gnuplot" << endl;
    }

#if defined( _MSC_VER )
    cin.get();
#endif

    return 0;
}
