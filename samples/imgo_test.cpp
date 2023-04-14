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

const int test_func_number = 6; // 0 - f1, 1 - f2, 2 - f3, ... 

double f1(double x, int j) {
    switch (j) {
        case 1: return exp(-sin(3.0 * x)) - 1.0 / 10.0 * pow(x - 1.0 / 2.0, 2) - 1.0;
        case 2: return -13.0 / 6.0 * x + sin(13.0 / 4.0 * (2.0 * x + 5.0)) - 53.0 / 12.0;
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f2(double x, int j) {
    switch (j) {
        case 1: return 1.0 / 20.0 - exp(-2.0 / 5.0 * (x + 5.0)) * sin(4.0 / 5.0 * M_PI * (x + 5.0));
        case 2: return (11.0 * x * x - 10.0 * x + 21.0) / (2.0 * (x * x + 1));
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f3(double x, int j) {
    double sum = 0.0;
    switch (j) {
        case 1: return 3.0 / 2.0 * (cos(7.0 / 20.0 * (x + 10.0)) - sin(7.0 / 4.0 * (x + 10.0)) + 1.0 / 2.0);
        case 2:
            for (int i = 1; i <= 5; i++) {
                sum += cos(i * x);
            }
            return -sum;
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f4(double x, int j) {
    double sum = 0.0;
    switch (j) {
        case 1:
            for (int i = 1; i <= 5; i++) {
                sum += cos(5.0 / 4.0 * (i + 1.0) * x + i);
            }
            return 6.0 / 25.0 - sum;
        case 2: return 9.0 / 50.0 - 9.0 / 2.0 * exp(-(x - 1.0 / 10.0)) * 
                sin(2.0 * M_PI * (x - 1.0 / 10.0));
        case 3: return 4.0 * sin(M_PI / 4.0 * x + 1.0 / 20.0) * 
                pow(pow(sin(M_PI / 2.0 * x + 1.0 / 10.0), 3.0) + pow(cos(M_PI / 2.0 * x + 1.0 / 10.0), 3.0), 2.0);
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f5(double x, int j) {
    switch (j) {
        case 1: return 17.0 / 25.0 - 2.0 / 29763.233 * (-1.0 / 6.0 * pow(x, 6) + 52.0 / 25.0 * pow(x, 5) - 39.0 / 80.0 * pow(x, 4) -
            71.0 / 10.0 * x * x * x + 79.0 / 20.0 * x * x + x - 1.0 / 10.0);
        case 2: return -14.0 / 125.0 * (3.0 * x - 8.0) * sin(252.0 / 125.0 * (x + 3.0 / 2.0)) - 1.0 / 2.0;
        case 3: return sin(0.423531 * x + 3.13531) + sin(10.0 / 3.0 * (0.423531 * x + 3.13531)) + log(0.423531 * x + 3.13531) + 0.36634 - 0.355766 * x;
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f6(double x, int j) {
    switch (j) {
        case 1: return 40.0 * (cos(4.0 * x) * (x - sin(x)) * exp( -(x * x) / 2.0));
        case 2: return 2.0 / 25.0 * (x + 4.0) - sin(12.0 / 5.0 * (x + 4.0));
        case 3: return -7.0 / 40.0 * (3.0 * x + 4.0) * sin(63.0 / 20.0 * (x + 4.0));
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f7(double x, int j) {
    switch (j) {
        case 1: return pow(sin(x), 3) * exp(-sin(3.0 * x)) + 1.0 / 2.0;
        case 2: return cos(7.0 / 5.0 * (x + 3.0)) - sin(7.0 * (x + 3.0)) + 3.0 / 10.0;
        case 3: return exp(-cos(4.0 * x - 3.0)) + 1.0 / 250.0 * (4.0 * x - 3.0) * (4.0 * x - 3.0) - 1.0;
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f8(double x, int j) {
    double sum = 0.0;
    switch (j) {
        case 1:
            return exp(-sin(4.0 * x)) - 1.0 / 10.0 * pow(x - 1.0 / 2.0, 2) - 1.0;
        case 2: 
            for (int i = 1; i <= 5; i++) {
                sum += cos(5.0 * (i + 1.0) * (x + 1.0 / 2.0));
            }
            return 3.0 / 10.0 - sum;
        case 3: return (-21.0 / 20.0 * x - 13.0 / 8.0) * sin(63.0 / 10.0 * x + 63.0 / 4.0) + 1.0 / 5.0;
        case 4: return cos(7.0 / 4.0 * x + 241.0 / 40.0) - sin(35.0 / 4.0 * x + 241.0 / 8.0) - 5.0;
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f9(double x, int j) {
    double sum = 0.0;
    switch (j) {
        case 1: return 1.0 / 40.0 * (x - 4.0) * (x - 32.0 / 5.0) * (x - 9.0) * (x - 11.0) *
            exp(-1.0 / 10.0 * pow(x - 13.0 / 2.0, 2));
        case 2: return (pow(sin(x + 1.0), 3) + pow(cos(x + 1.0), 3)) * exp(-(x + 1.0) / 10.0);
        case 3: return exp(-cos(3.0 / 5.0 * (x - 5.0 / 2.0))) + 1.0 / 10.0 * pow(3.0 / 25.0 * x - 4.0 / 5.0, 2) - 1.0;
        case 4:
            for (int i = 1; i <= 5; i++) {
                sum += 1.0 / 5.0 * sin((i + 1.0) * x - 1.0) + 2.0;
            }
            return sum;
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f10(double x, int j) {
    double a, b;
    switch (j) {
        case 1: return 2.0 * exp(-2.0 / M_PI * x) * sin(4.0 * x);
        case 2:
            a = 2.0 / M_PI * x - 1.0 / 2.0; 
            return -a * a * (-a * a + 5.0 * a - 6.0) / (a * a + 1.0) - 1.0 / 2.0;
        case 3: return pow(sin(x), 3) + pow(cos(2.0 * x), 3) - 3.0 / 10.0;
        case 4:
            b = 4.0 / M_PI * (x - 3.0 / 10.0) - 4.0;
            return -1.0 / 500.0 * b * b * b * b * b * b + 3.0 / 100.0 * b * b * b * b - 27.0 / 500.0 * b * b + 3.0 / 2.0;
        default: return numeric_limits<double>::quiet_NaN();
    }
}

int main() {
    ofstream ofstr("output_data/imgo_test.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstr_opt("output_data/imgo_test_opt.txt");
    if (!ofstr_opt.is_open()) cerr << "File opening error\n";

    double eps = 0.000001, r = 2.0, d = 0.0, x;
    int countIters, countEvals;
    int maxIters = 100000, maxEvals = 100000;

    vector<task_imgo> task_array = { task_imgo(f1, "f1(x)", 1, -2.5, 1.5, 1.05738,
                                     vector<double>{4.640837, 8.666667}, eps, maxIters, maxEvals, r, d),
                                     task_imgo(f2, "f2(x)", 1, -5.0, 5.0, 1.016,
                                     vector<double>{2.513269, 6.372595}, eps, maxIters, maxEvals, r, d),
                                     task_imgo(f3, "f3(x)", 1, -10.0, 10.0, -5.9921,
                                     vector<double>{3.124504, 13.201241}, eps, maxIters, maxEvals, r, d),
                                     task_imgo(f4, "f4(x)", 2, 0.0, 4.0, 2.45956,
                                     vector<double>{29.731102, 35.390605, 12.893183}, eps, maxIters, maxEvals, r, d),
                                     task_imgo(f5, "f5(x)", 2, -1.5, 11.0, 9.28491,
                                     vector<double>{5.654617, 0.931981, 2.021595}, eps, maxIters, maxEvals, r, d),
                                     task_imgo(f6, "f6(x)", 2, -4.0, 4.0, 2.32396,
                                     vector<double>{2.48, 25.108154, 8.835339}, eps, maxIters, maxEvals, r, d),
                                     task_imgo(f7, "f7(x)", 2, -3.0, 2.0, -0.774575,
                                     vector<double>{8.332010, 5.359309, 6.387862}, eps, maxIters, maxEvals, r, d),
                                     task_imgo(f8, "f8(x)", 3, -2.5, 1.5, -1.12724,
                                     vector<double>{20.184982, 90.598898, 6.372137, 10.415012}, eps, maxIters, maxEvals, r, d),
                                     task_imgo(f9, "f9(x)", 3, 0.0, 14.0, 4.0,
                                     vector<double>{0.873861, 1.682731, 1.254588, 3.843648}, eps, maxIters, maxEvals, r, d),
                                     task_imgo(f10, "f10(x)", 3, 0.0, 2.0 * M_PI, 4.2250023,
                                     vector<double>{3.170468, 4.329008, 7.999984, 12.442132}, eps, maxIters, maxEvals, r, d) };

    imgo_method imgo(nullptr);

    vector<double> lambdas;
    vector<trial_constr> trials;
    for (int i = 0; i < task_array.size(); i++) {
        if (task_array[i].used) {
            imgo.setF(task_array[i].f);
            imgo.setNumberConstraints(task_array[i].numberConstraints);
            imgo.setAB(task_array[i].A[0], task_array[i].B[0]);
            imgo.setEps(task_array[i].eps);
            imgo.setMaxIters(task_array[i].maxIters);
            imgo.setMaxEvals(task_array[i].maxEvals);
            imgo.setR(task_array[i].r);
            imgo.setD(task_array[i].d);

            imgo.solve(countIters, countEvals, x);
            imgo.getLambda(lambdas);

            cout << "Function: " << task_array[i].name << endl;
            cout << "Number of constraints = " << task_array[i].numberConstraints << endl;
            cout << "[a; b] = [" << task_array[i].A[0] << "; " << task_array[i].B[0] << "]"<< endl;
            cout << "Lipschitz constant:" << endl;
            cout << "L*(" << task_array[i].name << ") = " << task_array[i].L[task_array[i].numberConstraints] << endl;
            for (int j = 0; j < task_array[i].numberConstraints; j++) {
                cout << "L*(g" << j + 1 << ") = " << task_array[i].L[j] << endl;
            }
            cout << "X* = " << setprecision(8) << task_array[i].X_opt[0] << endl;
            cout << "f(X*) = " << setprecision(8) << task_array[i].f(task_array[i].X_opt[0], 
                                                                     task_array[i].numberConstraints + 1) << endl;
            cout << "Parameters for method:" << endl;
            cout << "eps = " << eps << " r = " << r << " d = " << d << endl;
            cout << "Trials result:" << endl;
            cout << "Number of trials = " << countIters << endl;
            cout << "Number of evals = " << countEvals << endl;
            cout << "Estimation of the Lipschitz constant:" << endl;
            cout << "L(" << task_array[i].name << ") = " << lambdas[task_array[i].numberConstraints] << endl;
            for (int j = 0; j < task_array[i].numberConstraints; j++) {
                cout << "L(g" << j + 1 << ") = " << lambdas[j] << endl;
            }
            cout << "X = " << setprecision(8) << x << endl;
            cout << "f(X) = " << setprecision(8) << task_array[i].f(x, task_array[i].numberConstraints + 1) << endl;
            cout << "|X* - X| = " << setprecision(8) << abs(task_array[i].X_opt[0] - x) << endl;
            cout << "|f(X*) - f(X)| = " << setprecision(8) << abs(task_array[i].f(task_array[i].X_opt[0], 
                                                                  task_array[i].numberConstraints + 1) - 
                                                                  task_array[i].f(x, task_array[i].numberConstraints + 1))
                                                                                                                   << endl;
            cout << endl;

            // Saving points for plotting
            ofstr << task_array[i].X_opt[0] << " " << task_array[i].f(task_array[i].X_opt[0],
                                                                      task_array[i].numberConstraints + 1) << endl;
            ofstr << endl << endl;
            ofstr << x << " " << task_array[i].f(x, task_array[i].numberConstraints + 1) << endl;
            ofstr << endl << endl;
            imgo.getTrialPoints(trials);
            for (int j = 0; j < trials.size(); j++) {
                ofstr << trials[j].x << " " << trials[j].z << endl;
            }
            ofstr << endl << endl;
        }
    }
    ofstr.close();

    ofstr_opt << "array functions[" << task_array.size() << "]" << endl;
    for (int i = 0; i < task_array.size(); i++) {
        ofstr_opt << "functions[" << i + 1 << "] = \"f" << i + 1 << "(x) title \\\"f(x)\\\"";
        for (int j = 0; j < task_array[i].numberConstraints; j++) {
            ofstr_opt << ", g" << i + 1 << "_" << j + 1 << "(x) title \\\"g" << j + 1 <<"(x)\\\"";
        }
        ofstr_opt << "\""<< endl;
    }
    ofstr_opt.close();

    // Plotting the function(works with gnuplot)
    int error;
#if defined(__linux__)
    setenv("QT_QPA_PLATFORM", "xcb", false);
    error = system("chmod +x scripts/imgo_test.gp");
    if (error != 0) {
        cerr << "Error chmod" << endl;
    }
#endif

    char str[100];
    sprintf(str, "gnuplot -c scripts/imgo_test.gp %d", test_func_number);
    error = system(str);
    if (error != 0) {
        cerr << "Error gnuplot" << endl;
    }

#if defined( _MSC_VER )
    cin.get();
#endif

    return 0;
}