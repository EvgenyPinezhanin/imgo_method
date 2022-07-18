#include <iostream>
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
    return -(-sin(x) -sin(10.0 / 3.0 * x) - log(x) + 0.84 * x - 3);
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
    return -((-x * x + 5.0 * x -6.0) / (x * x + 1));
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
    return -(-2.0 * log(x - 2.0) - 1);
}

double f19(double x) {
    return -(x - sin(3.0 * x) + 1);
}

double f20(double x) {
    return -(x - sin(x)) * exp(- x * x);
}

int main() {
    double x, eps = 0.001, r = 2.0; // > 1
    int count, Nmax = 1000;
    Stop stop = Stop::ACCURACY;

    vector<task_gsa> task_array = { task_gsa(f1, "f1(x)", -1.5, 11.0, 10.0, eps, Nmax, r, stop),
                                    task_gsa(f2, "f2(x)", 2.7, 7.5, 5.145735, eps, Nmax, r, stop),
                                    task_gsa(f3, "f3(x)", -10.0, 10.0, 5.791785, eps, Nmax, r, stop),
                                    task_gsa(f4, "f4(x)", 1.9, 3.9, 2.868, eps, Nmax, r, stop),
                                    task_gsa(f5, "f5(x)", 0.0, 1.2, 0.96609, eps, Nmax, r, stop),
                                    task_gsa(f6, "f6(x)", -10.0, 10.0, 0.67956, eps, Nmax, r, stop),
                                    task_gsa(f7, "f7(x)", 2.7, 7.5, 5.19978, eps, Nmax, r, stop),
                                    task_gsa(f8, "f8(x)", -10.0, 10.0, -7.0835, eps, Nmax, r, stop),
                                    task_gsa(f9, "f9(x)", 3.1, 20.4, 17.039, eps, Nmax, r, stop),
                                    task_gsa(f10, "f10(x)", 0.0, 10.0, 7.9787, eps, Nmax, r, stop),
                                    task_gsa(f11, "f11(x)", -1.57, 6.28, 2.094, eps, Nmax, r, stop),
                                    task_gsa(f12, "f12(x)", 0.0, 6.28, 3.142, eps, Nmax, r, stop),
                                    task_gsa(f13, "f13(x)", 0.001, 0.99, 0.7071, eps, Nmax, r, stop),
                                    task_gsa(f14, "f14(x)", 0.0, 4.0, 0.224885, eps, Nmax, r, stop),
                                    task_gsa(f15, "f15(x)", -5.0, 5.0, 2.4142, eps, Nmax, r, stop),
                                    task_gsa(f16, "f16(x)", -3.0, 3.0, 3.0, eps, Nmax, r, stop),
                                    task_gsa(f17, "f17(x)", -4.0, 4.0, -3.0, eps, Nmax, r, stop),
                                    task_gsa(f18, "f18(x)", 0.0, 6.0, 2.0, eps, Nmax, r, stop),
                                    task_gsa(f19, "f19(x)", 0.0, 6.5, 5.87287, eps, Nmax, r, stop),
                                    task_gsa(f20, "f20(x)", -10.0, 10.0, 1.195137, eps, Nmax, r, stop) };

    gsa_method gsa(nullptr);

    for (int i = 0; i < task_array.size(); i++) {
        if (task_array[i].used) {
            gsa.setF(task_array[i].f);
            gsa.setAB(task_array[i].A[0], task_array[i].B[0]);
            gsa.setEps(task_array[i].eps);
            gsa.setNmax(task_array[i].Nmax);
            gsa.setR(task_array[i].r);

            gsa.solve(count, x, task_array[i].stop);

            cout << "Function: " << task_array[i].name << endl;
            cout << "[a; b] = [" << task_array[i].A[0] << "; " << task_array[i].B[0] << "]"<< endl;
            cout << "X* = " << setprecision(8) << task_array[i].X_opt[0] << endl;
            cout << "X = " << setprecision(8)  << x << endl;
            cout << "|X* - X| = " << setprecision(8) << abs(task_array[i].X_opt[0] - x) << endl;
            cout << "|f(X*) - f(X)| = " << setprecision(8) << abs(task_array[i].f(task_array[i].X_opt[0]) - task_array[i].f(x)) << endl;
            cout << "Count of trials = " << count << endl;
            cout << endl;
        }
    }

#if defined( _MSC_VER )
    cin.get();
#endif
    return 0;
}
