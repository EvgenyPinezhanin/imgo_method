/* #include <iostream>
#include <fstream> 
#include <cmath>
#include <Grishagin/grishagin_function.hpp>
#include <GKLS/GKLSProblem.hpp>
#include <imgo.h>
#include <task.h>

using namespace std;

const int test_func_number = 1; // 0, 1 

double f_test_1(vector<double> x, int j) {
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

double f_test_2(vector<double> x, int j) {
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
    ofstream ofstr("peano_test_trial_points.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";

    int n = 2, m = 0, den = 10, key = 1, Nmax = 5000;
    double eps = 0.001, r = 2.0, d = 0.001;
    Stop stop = ACCURACY;
    vector<double> X(2);
    vector<vector<double>> trial_vec;
    int number_trials;

    vector<Task_peano> task{ Task_peano(f_test_1, "Test-1", n, 3, vector<double>{0.0, -1.0}, vector<double>{4.0, 3.0},
                                        vector<double>{0.942, 0.944}, eps, Nmax, r, d, den, key, stop, 1),
                             Task_peano(f_test_2, "Test-2", n, 4, vector<double>{0.0, 0.0}, vector<double>{80.0, 80.0},
                                        vector<double>{77.489, 63.858}, 0.001, Nmax, 3.3, 0.01, den, key, stop, 1) };

    imgo_method imgo(nullptr, 2, 0, vector<double>{0.0, 0.0}, vector<double>{1.0, 1.0});

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

            // Подготовка данных для построения графика
            ofstr << X[0] << " " << X[1] << " " << task[i].f(X, task[i].m + 1) << endl;
            ofstr << endl << endl;
            ofstr << task[i].X_opt[0] << " " << task[i].X_opt[1] << " " << task[i].f(task[i].X_opt, task[i].m + 1) << endl;
            ofstr << endl << endl;
            imgo.getPoints(trial_vec);
            for (int j = 0; j < trial_vec.size(); j++) {
                ofstr << trial_vec[j][0] << " " << trial_vec[j][1] << " " << task[i].f(trial_vec[j], task[i].m + 1) << endl;
            }
            ofstr << endl << endl;
        }
    }
    ofstr.close();

    // Построение графика(работает с помощью gnuplot)
#if defined(__linux__)
    int error;
    setenv("QT_QPA_PLATFORM", "xcb", false);
    error = system("chmod +x scripts/peano_test.gp");
    if (error != 0) {
        cerr << "Error chmod" << std::endl;
    }

    char str[100];
    sprintf(str, "gnuplot -p -c scripts/peano_test.gp %d", test_func_number);
    error = system(str);
    if (error != 0) {
        cerr << "Error gnuplot" << std::endl;
    }
#endif

    return 0;
}
 */