/* #include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <imgo.h>
#include <task.h>

using namespace std;

const int sample_func_number = 3; // 0, 1, 2, 3

double f1_sample(vector<double> x, int j) {
    switch (j) {
        case 1: return 1.0 - x[0] - x[1];
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f2_sample(vector<double> x, int j) {
    switch (j) {
        case 1: return (x[0] - 1.0) * (x[0] - 1.0) / 5.0 + (x[1] - 1.0) * (x[1] - 1.0) / 5.0;
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f3_sample(vector<double> x, int j) {
    switch (j) {
        case 1: return 1.0 - x[0] - x[1];
        case 2: return x[0] * x[0] / 5.0 + x[1] * x[1] / 5.0;
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f4_sample(vector<double> x, int j) {
    switch (j) {
        case 1: return pow(x[0] - 2.0, 2) + pow(x[1] - 2.0, 2) - 2.0;
        case 2: return x[0] * x[0] / 5.0 + x[1] * x[1] / 5.0;
        default: return numeric_limits<double>::quiet_NaN();
    }
}

int main() {
    ofstream ofstr("peano_sample_trial_points.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";

    vector<double> X(2);
    std::vector<vector<double>> trial_vec;
    int number_trials;
    int n = 2, den = 10, key = 1, Nmax = 1000;
    double eps = 0.001, r = 2.0, d = 0.0;
    Stop stop = ACCURACY;

    vector<Task_peano> task = { Task_peano(f1_sample, "f1_sample", n, 0, vector<double>{-4.0, -4.0}, vector<double>{4.0, 4.0},
                                           vector<double>{4.0, 4.0}, eps, Nmax, r, d, den, key, stop, 1),
                                Task_peano(f2_sample, "f2_sample", n, 0, vector<double>{-4.0, -4.0}, vector<double>{4.0, 4.0},
                                           vector<double>{1.0, 1.0}, eps, Nmax, r, d, den, key, stop, 1),
                                Task_peano(f3_sample, "f3_sample", n, 1, vector<double>{-1.0, -1.0}, vector<double>{1.0, 1.0},
                                           vector<double>{0.5, 0.5}, eps, Nmax, r, d, den, key, stop, 1),
                                Task_peano(f4_sample, "f4_sample", n, 1, vector<double>{0.0, 0.0}, vector<double>{3.0, 3.0},
                                           vector<double>{1.0, 1.0}, eps, Nmax, r, d, den, key, stop, 1) };

    imgo_method imgo(nullptr, n, 0, vector<double>{0.0, 0.0}, vector<double>{0.0, 0.0}, r, d, eps, Nmax, den, key);

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
            cout << "n = " << task[i].n << " m = " << task[i].den << " key = " << task[i].key << endl;
            cout << "Trials result:" << endl;
            cout << "Number of trials = " << number_trials << endl;
            cout << "x*_min = " << task[i].X_opt[0] << " y*_min = " << task[i].X_opt[1] << endl;
            cout << "x_min = " << X[0] << " y_min = " << X[1] << endl;
            cout << "||X* - X|| = " << sqrt((task[i].X_opt[0] - X[0]) * (task[i].X_opt[0] - X[0]) + 
                                            (task[i].X_opt[1] - X[1]) * (task[i].X_opt[1] - X[1])) << endl;
            cout << "|f(X*) - f(X)| = " << abs(task[i].f(task[i].X_opt, task[i].m + 1) - task[i].f(X, task[i].m + 1)) << endl;
            cout << endl;

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

    // Построение графика(работает только под Lunux с помощью gnuplot)
#if defined(__linux__)
    int error;
    setenv("QT_QPA_PLATFORM", "xcb", false);
    error = system("chmod +x scripts/peano_sample.gp");
    if (error != 0) {
        std::cerr << "Error chmod" << std::endl;
    }

    char str[100];
    sprintf(str, "gnuplot -p -c scripts/peano_sample.gp %d", sample_func_number);
    error = system(str);
    if (error != 0) {
        std::cerr << "Error gnuplot" << std::endl;
    }
#endif

    return 0;
}
 */