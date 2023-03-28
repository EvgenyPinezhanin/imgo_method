#if defined( _MSC_VER )
    #define _CRT_SECURE_NO_WARNINGS
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#if defined( _MSC_VER )
    #define _USE_MATH_DEFINES
    #include <math.h>
#else
    #include <cmath>
#endif

#include <gsa.h>
#include <task.h>

using namespace std;

const int test_func_number = 15; // 0 - f1, 1 - f2, 2 - f3, 3 - f4, ...

double f1(double x) {
    return -(-1.0 / 6.0 * pow(x, 6) + 52.0 / 25.0 * pow(x, 5) - 39.0 / 80.0 * pow(x, 4) - 
             71.0 / 10.0 * pow(x, 3) + 79.0 / 20.0 * x * x + x - 1.0 / 10.0);
}

double f2(double x) {
    return -(-sin(x) - sin(10.0 / 3.0 * x));
}

double f3(double x) {
    double sum = 0;
    for (int i = 1; i <= 5; i++) {
        sum += i * sin((i + 1.0) * x  + i);
    }
    return -sum;
}

double f4(double x) {
    return -(16.0 * x * x - 24.0 * x + 5.0) * exp(-x);
}

double f5(double x) {
    return -(-3.0 * x + 1.4) * sin(18.0 * x);
}

double f6(double x) {
    return -((x + sin(x)) * exp(-x * x));
}

double f7(double x) {
    return -(-sin(x) - sin(10.0 / 3.0 * x) - log(x) + 0.84 * x - 3.0);
}

double f8(double x) {
    double sum = 0;
    for (int i = 1; i <= 5; i++) {
        sum += i * cos((i + 1.0) * x  + i);
    }
    return -sum;
}

double f9(double x) {
    return -(-sin(x) - sin(2.0 / 3.0 * x));
}

double f10(double x) {
    return -(x * sin(x));
}

double f11(double x) {
    return -(-2.0 * cos(x) - cos(2.0 * x));
}

double f12(double x) {
    return -(-pow(sin(x), 3) - pow(cos(x), 3));
}

double f13(double x) {
    if (x * x - 1 < 0) {
        return -(pow(x, 2.0 / 3.0) + pow(-(x * x - 1), 1.0 / 3.0));
    }
    return -(pow(x, 2.0 / 3.0) - pow(x * x - 1, 1.0 / 3.0));
}

double f14(double x) {
    return -(exp(-x) * sin(2 * M_PI * x));
}

double f15(double x) {
    return -((-x * x + 5.0 * x - 6.0) / (x * x + 1));
}

double f16(double x) {
    return -(-2.0 * (x - 3) * (x - 3) - exp(- x * x / 2));
}

double f17(double x) {
    return -(-pow(x, 6) + 15.0 * pow(x, 4) - 27.0 * x * x - 250.0);
}

double f18(double x) {
    if (x <= 3.0) {
        return (x - 2.0) * (x - 2.0);
    } 
    return -(-2.0 * log(x - 2.0) - 1.0);
}

double f19(double x) {
    return -(x - sin(3.0 * x) + 1.0);
}

double f20(double x) {
    return -(x - sin(x)) * exp(- x * x);
}

int main() {
    ofstream ofstr("output_data/gsa_test.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";

    double x, eps = 0.001, r = 2.0; // > 1
    int countIters, countTrials, countEvals;
    int maxIters = 100000, maxEvals = 100000;

    vector<task_gsa> task_array = { task_gsa(f1, "f1(x)", -1.5, 11.0, 10.0, 13870.0, eps, maxIters, maxEvals, r),
                                    task_gsa(f2, "f2(x)", 2.7, 7.5, 5.145735, 4.29, eps, maxIters, maxEvals, r),
                                    task_gsa(f3, "f3(x)", -10.0, 10.0, -0.49139, 67.0, eps, maxIters, maxEvals, r),
                                    task_gsa(f4, "f4(x)", 1.9, 3.9, 2.868, 3.0, eps, maxIters, maxEvals, r),
                                    task_gsa(f5, "f5(x)", 0.0, 1.2, 0.96609, 36.0, eps, maxIters, maxEvals, r),
                                    task_gsa(f6, "f6(x)", -10.0, 10.0, 0.67956, 2.5, eps, maxIters, maxEvals, r),
                                    task_gsa(f7, "f7(x)", 2.7, 7.5, 5.19978, 6.0, eps, maxIters, maxEvals, r),
                                    task_gsa(f8, "f8(x)", -10.0, 10.0, 5.48286, 67.0, eps, maxIters, maxEvals, r),
                                    task_gsa(f9, "f9(x)", 3.1, 20.4, 17.039, 1.7, eps, maxIters, maxEvals, r),
                                    task_gsa(f10, "f10(x)", 0.0, 10.0, 7.9787, 11.0, eps, maxIters, maxEvals, r),
                                    task_gsa(f11, "f11(x)", -1.57, 6.28, 2.094, 3.0, eps, maxIters, maxEvals, r),
                                    task_gsa(f12, "f12(x)", 0.0, 6.28, 4.712, 2.2, eps, maxIters, maxEvals, r),
                                    task_gsa(f13, "f13(x)", 0.001, 0.99, 0.7071, 8.5, eps, maxIters, maxEvals, r),
                                    task_gsa(f14, "f14(x)", 0.0, 4.0, 0.224885, 6.5, eps, maxIters, maxEvals, r),
                                    task_gsa(f15, "f15(x)", -5.0, 5.0, 2.4142, 6.5, eps, maxIters, maxEvals, r),
                                    task_gsa(f16, "f16(x)", -3.0, 3.0, 3.0, 85.0, eps, maxIters, maxEvals, r),
                                    task_gsa(f17, "f17(x)", -4.0, 4.0, -3.0, 2520.0, eps, maxIters, maxEvals, r),
                                    task_gsa(f18, "f18(x)", 0.0, 6.0, 2.0, 4.0, eps, maxIters, maxEvals, r),
                                    task_gsa(f19, "f19(x)", 0.0, 6.5, 5.87287, 4.0, eps, maxIters, maxEvals, r),
                                    task_gsa(f20, "f20(x)", -10.0, 10.0, 1.195137, 1.3, eps, maxIters, maxEvals, r) };

    gsa_method gsa(nullptr);

    vector<trial> trials;
    for (int i = 0; i < task_array.size(); i++) {
        if (task_array[i].used) {
            gsa.setF(task_array[i].f);
            gsa.setAB(task_array[i].A[0], task_array[i].B[0]);
            gsa.setEps(task_array[i].eps);
            gsa.setMaxIters(task_array[i].maxIters);
            gsa.setMaxEvals(task_array[i].maxEvals);
            gsa.setR(task_array[i].r);

            gsa.solve(countIters, countTrials, countEvals, x);

            cout << "Function: " << task_array[i].name << endl;
            cout << "[a; b] = [" << task_array[i].A[0] << "; " << task_array[i].B[0] << "]"<< endl;
            cout << "Lipschitz constant = " << task_array[i].L[0] << endl;
            cout << "X* = " << setprecision(8) << task_array[i].X_opt[0] << endl;
            cout << "f(X*) = " << setprecision(8) << task_array[i].f(task_array[i].X_opt[0]) << endl;
            cout << "Parameters for method:" << endl;
            cout << "eps = " << eps << " r = " << r << endl;
            cout << "Trials result:" << endl;
            cout << "Number of iters = " << countIters << endl;
            cout << "Number of trials = " << countTrials << endl;
            cout << "Number of evals = " << countEvals << endl;
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
    error = system("chmod +x scripts/gsa_test.gp");
    if (error != 0) {
        cerr << "Error chmod" << endl;
    }
#endif

    char str[100];
    sprintf(str, "gnuplot -c scripts/gsa_test.gp %d", test_func_number);
    error = system(str);
    if (error != 0) {
        cerr << "Error gnuplot" << endl;
    }

#if defined( _MSC_VER )
    cin.get();
#endif

    return 0;
}
