#include <iostream>
#include <iomanip>
#include <fstream>
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
    switch (j) {
        case 1: return exp(-sin(3.0 * x)) - 1.0 / 10.0 * pow(x - 1.0 / 2.0, 2.0) - 1.0;
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

// Error???
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

int main() {
    ofstream ofstr("output_data_test.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    vector<trial> trial_vec;

    double eps = 0.001, r = 3.0, d = 0.0, x_min;
    int count, Nmax = 1000;
    Stop stop = ACCURACY;

    vector<Task> task_array = { Task(f1, "f1(x)", 1, -2.5, 1.5, 1.05738, eps, Nmax, r, d, stop, 1),
                                Task(f2, "f2(x)", 1, -5.0, 5.0, 1.016, eps, Nmax, r, d, stop, 1),
                                Task(f4, "f4(x)", 2, 0.0, 4.0, 2.45956, eps ,Nmax, r, d, stop,  1),
                                Task(f5, "f5(x)", 2, -1.5, 11.0, 8.85725, eps, Nmax, r, d, stop, 1),
                                Task(f6, "f6(x)", 2, -4.0, 4.0, 2.32396, eps, Nmax, r, d, stop, 1),
                                Task(f8, "f8(x)", 3, -2.5, 1.5, -1.12724, eps, Nmax, r, d, stop, 1),
                                Task(f9, "f9(x)", 3, 0.0, 14.0, 4.0, eps, Nmax, r, d, stop, 1) };

    imgo_method imgo(nullptr, 0, 0.0, 0.0, r, d, eps);

    for (int i = 0; i < task_array.size(); i++) {
        imgo.setFunc(task_array[i].f);
        imgo.setM(task_array[i].m);
        imgo.setAB(task_array[i].A, task_array[i].B);
        imgo.setEps(task_array[i].eps);
        imgo.setNmax(task_array[i].Nmax);
        imgo.setR(task_array[i].r);
        imgo.setD(task_array[i].d);
        imgo.setNmax(task_array[i].Nmax);

        x_min = imgo.solve(count, task_array[i].stop);

        cout << "Function: " << task_array[i].name << endl;
        cout << "[a; b] = [" << task_array[i].A[0] << "; " << task_array[i].B[0] << "]"<< endl;
        cout << "X* = " << setprecision(12) << task_array[i].X_opt[0] << endl;
        cout << "X = " << setprecision(12) << x_min << endl;
        cout << "|X* - X| = " << abs(task_array[i].X_opt[0] - x_min) << endl;
        cout << "|f(X*) - f(X)| = " << abs(task_array[i].f(task_array[i].X_opt[0], task_array[i].m + 1) - 
                                           task_array[i].f(x_min, task_array[i].m + 1)) << endl;
        cout << "Number of trials = " << count << endl;
        cout << endl;

        imgo.getTrialPoints(trial_vec);
        ofstr << task_array[i].A[0] << " " << task_array[i].B[0] << " " << task_array[i].m 
              << " " << x_min << " " << task_array[i].X_opt[0] << endl;
        for (int j = 0; j < trial_vec.size(); j++) {
            ofstr << trial_vec[j].x << " " << trial_vec[j].z << endl;
        }
        ofstr << endl;
    }

    ofstr.close();
#if defined( _MSC_VER )
    cin.get();
#endif
    return 0;
}