#if defined( _MSC_VER )
    #define _CRT_SECURE_NO_WARNINGS    
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include <mggsa.h>
#include <task.h>

using namespace std;

const int test_func_number = 0; // 0 - f1, 1 - f2

double f1(vector<double> x, int j) {
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

const double C[20] = {75.1963666677,-3.8112755343,0.1269366345,-0.0020567665,0.000010345,
                      -6.8306567631,0.0302344793,-0.0012813448,0.0000352559,-0.0000002266,
                      0.2564581253,-0.0034604030,0.0000135139,-28.1064434908,-0.0000052375,
                      -0.0000000063,0.0000000007,0.0003405462,-0.0000016638,-2.8673112392 };

double f2(vector<double> x, int j) {
    switch (j) {
        case 1: return 450.0 - x[0] * x[1];
        case 2: return (0.1 * x[0] - 1.0) * (0.1 * x[0] - 1.0) - x[1];
        case 3: return 8.0 * (x[0] - 40.0) - (x[1] - 30.0) * (x[1] - 55.0);
        case 4: return x[1] + (x[0] - 35.0) * (x[0] - 30.0) / 125.0 - 80.0;
        case 5: return -(C[0] + C[1] * x[0] + C[2] * x[0] * x[0] + C[3] * pow(x[0], 3) + C[4] * pow(x[0], 4) + C[5] * x[1] +
                         C[6] * x[0] * x[1] + C[7] * x[0] * x[0] * x[1] + C[8] * pow(x[0], 3) * x[1] + C[9] * pow(x[0], 4) * x[1] +
                         C[10] * x[1] * x[1] + C[11] * pow(x[1], 3) + C[12] * pow(x[1], 4) + C[13] / (x[1] + 1) + C[14] * x[0] * x[0] * x[1] * x[1] +
                         C[15] * pow(x[0], 3) * x[1] * x[1] + C[16] * pow(x[0], 3) * pow(x[1], 3) + C[17] * x[0] * x[1] * x[1] +
                         C[18] * x[0] * pow(x[1], 3) + C[19] * exp(0.0005 * x[0] * x[1]));
        default: return numeric_limits<double>::quiet_NaN();
    }
}

int main() {
    ofstream ofstr("output_data/mggsa_test.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstr_opt("output_data/mggsa_test_opt.txt");
    if (!ofstr_opt.is_open()) cerr << "File opening error\n";

    double eps = 0.001, r = 2.2, d = 0.05;
    int count_trials, n = 2, den = 10, key = 1, Nmax = 1000;
    vector<double> X(n);
    Stop stop = Stop::ACCURACY;

    vector<task_mggsa> task_array = { task_mggsa(f1, "f1", n, 3, vector<double>{0.0, -1.0}, vector<double>{4.0, 3.0},
                                                 vector<double>{0.942, 0.944}, vector<double>{}, eps, Nmax, r, d, 12, key, stop, true),
                                      task_mggsa(f2, "f2", n, 4, vector<double>{0.0, 0.0}, vector<double>{80.0, 80.0},
                                                 vector<double>{77.489, 63.858}, vector<double>{}, eps, Nmax, 3.3, 0.01, den, key, stop, true) };

    mggsa_method mggsa(nullptr);

    vector<double> mu;
    vector<vector<double>> points;
    for (int i = 0; i < task_array.size(); i++) {
        if (task_array[i].used) {
            mggsa.setF(task_array[i].f);
            mggsa.setN(task_array[i].n);
            mggsa.setM(task_array[i].m);
            mggsa.setAB(task_array[i].A, task_array[i].B);
            mggsa.setEps(task_array[i].eps);
            mggsa.setNmax(task_array[i].Nmax);
            mggsa.setR(task_array[i].r);
            mggsa.setD(task_array[i].d);
            mggsa.setDen(task_array[i].den);
            mggsa.setKey(task_array[i].key);

            mggsa.solve(count_trials, X, task_array[i].stop);
            mggsa.getLambda(mu);

            cout << "Function: " << task_array[i].name << endl;
            cout << "Dimension = " << task_array[i].n << endl;
            cout << "Number of constrained = " << task_array[i].m << endl;
            cout << "[A; B] = [(" << task_array[i].A[0] << ", " << task_array[i].A[1] << "); (" << 
                                     task_array[i].B[0] << ", " << task_array[i].B[1] << ")]"<< endl;
            cout << "X* = (" << task_array[i].X_opt[0] << ", " << task_array[i].X_opt[1] << ")" << endl;
            cout << "f(X*) = " << task_array[i].f(task_array[i].X_opt, task_array[i].m + 1) << endl;
            cout << "Parameters for method:" << endl;
            cout << "eps = " << task_array[i].eps << " r = " << task_array[i].r << 
                    " d = " << task_array[i].d << endl;
            cout << "Parameters for constructing the Peano curve:" << endl;
            cout << "m = " << task_array[i].den << " key = " << task_array[i].key << endl;
            cout << "Trials result:" << endl;
            cout << "Number of trials = " << count_trials << endl;
            cout << "Estimation of the Lipschitz constant:" << endl;
            cout << "L(" << task_array[i].name << ") = " << mu[task_array[i].m] << endl;
            for (int j = 0; j < task_array[i].m; j++) {
                cout << "L(g" << j + 1 << ") = " << mu[j] << endl;
            }
            cout << "X = (" << X[0] << ", " << X[1] << ")" << endl;
            cout << "f(X) = " << task_array[i].f(X, task_array[i].m + 1) << endl;
            cout << "||X* - X|| = " << sqrt((task_array[i].X_opt[0] - X[0]) * (task_array[i].X_opt[0] - X[0]) + 
                                            (task_array[i].X_opt[1] - X[1]) * (task_array[i].X_opt[1] - X[1])) << endl;
            cout << "|f(X*) - f(X)| = " << abs(task_array[i].f(task_array[i].X_opt, task_array[i].m + 1) - 
                                               task_array[i].f(X, task_array[i].m + 1)) << endl;
            cout << endl;

            // Saving points for plotting
            ofstr << task_array[i].X_opt[0] << " " << task_array[i].X_opt[1] << " " << 
                     task_array[i].f(task_array[i].X_opt, task_array[i].m + 1) << endl;
            ofstr << endl << endl;
            ofstr << X[0] << " " << X[1] << " " << task_array[i].f(X, task_array[i].m + 1) << endl;
            ofstr << endl << endl;
            mggsa.getPoints(points);
            for (int j = 0; j < points.size(); j++) {
                ofstr << points[j][0] << " " << points[j][1] << " " << 
                         task_array[i].f(points[j], task_array[i].m + 1) << endl;
            }
            ofstr << endl << endl;
        }
    }
    ofstr.close();

    size_t size = task_array.size();
    ofstr_opt << "array AX[" << size << "]" << endl;
    ofstr_opt << "array AY[" << size << "]" << endl;
    ofstr_opt << "array BX[" << size << "]" << endl;
    ofstr_opt << "array BY[" << size << "]" << endl;
    for (int i = 0; i < size; i++) {
        ofstr_opt << "AX[" << i + 1 << "] = " << task_array[i].A[0] << endl;
        ofstr_opt << "BX[" << i + 1 << "] = " << task_array[i].B[0] << endl;
        ofstr_opt << "AY[" << i + 1 << "] = " << task_array[i].A[1] << endl;
        ofstr_opt << "BY[" << i + 1 << "] = " << task_array[i].B[1] << endl;
    }
    ofstr_opt.close();

    // Plotting the function(works with gnuplot)
    int error;
#if defined(__linux__)
    setenv("QT_QPA_PLATFORM", "xcb", false);
    error = system("chmod +x scripts/mggsa_test.gp");
    if (error != 0) {
        cerr << "Error chmod" << endl;
    }
#endif

    char str[100];
    sprintf(str, "gnuplot -c scripts/mggsa_test.gp %d", test_func_number);
    error = system(str);
    if (error != 0) {
        cerr << "Error gnuplot" << endl;
    }

#if defined(_MSC_VER)
    cin.get();
#endif

    return 0;
}
