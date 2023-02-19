#include <iostream>
#include <vector>
#include <cmath>

#include <direct_method.h>
#include <task.h>

using namespace std;

double f1(int n, const double *X, int *undefined_flag, void *data) {
    int *count = static_cast<int*>(data); (*count)++;
    return 1.0 - X[0] - X[1];
}

double f2(int n, const double *X, int *undefined_flag, void *data) {
    int *count = static_cast<int*>(data); (*count)++;
    return (X[0] - 1.0) * (X[0] - 1.0) / 5.0 + (X[1] - 1.0) * (X[1] - 1.0) / 5.0;
}

double f3(int n, const double *X, int *undefined_flag, void *data) {
    int *count = static_cast<int*>(data); (*count)++;
    if (1.0 - X[0] - X[1] > 0.0) *undefined_flag = 1;
    return X[0] * X[0] / 5.0 + X[1] * X[1] / 5.0;
}

double f4(int n, const double *X, int *undefined_flag, void *data) {
    int *count = static_cast<int*>(data); (*count)++;
    if ((X[0] - 2.0) * (X[0] - 2.0) + (X[1] - 2.0) * (X[1] - 2.0) - 2.0 > 0.0) *undefined_flag = 1;
    return X[0] * X[0] / 5.0 + X[1] * X[1] / 5.0;
}

int main() {
    int count_trials, n = 2;
    int max_feval = 10000, max_iter = 10000;
    double magic_eps = 1.0e-4, eps = 0.01;
    double volume_reltol = eps / sqrt(n);
    double sigma_reltol  = 0.0;
    direct_algorithm algorithm = DIRECT_ORIGINAL;
    vector<double> X;
    double minf;

    vector<task_direct> task_array = { task_direct("f1", f1, &count_trials, n, vector<double>{-4.0, -4.0}, vector<double>{4.0, 4.0},
                                                   vector<double>{4.0, 4.0}, vector<double>{}, max_feval, max_iter, magic_eps,
                                                   volume_reltol, sigma_reltol, nullptr, algorithm),
                                       task_direct("f2", f2, &count_trials, n, vector<double>{-4.0, -4.0}, vector<double>{4.0, 4.0},
                                                   vector<double>{1.0, 1.0}, vector<double>{}, max_feval, max_iter, magic_eps,
                                                   volume_reltol, sigma_reltol, nullptr, algorithm),
                                       task_direct("f3", f3, &count_trials, n, vector<double>{-1.0, -1.0}, vector<double>{1.0, 1.0},
                                                   vector<double>{0.5, 0.5}, vector<double>{}, max_feval, max_iter, magic_eps,
                                                   volume_reltol, sigma_reltol, nullptr, algorithm),
                                       task_direct("f4", f4, &count_trials, n, vector<double>{0.0, 0.0}, vector<double>{3.0, 3.0},
                                                   vector<double>{1.0, 1.0}, vector<double>{}, max_feval, max_iter, magic_eps,
                                                   volume_reltol, sigma_reltol, nullptr, algorithm) };

    direct_method direct;

    int flag_tmp, f_data_tmp;
    int *count_tmp;
    for (int i = 0; i < task_array.size(); i++) {
        if (task_array[i].used) {
            count_tmp = (static_cast<int*>(task_array[i].f_data));
            *count_tmp = 0;

            direct.setF(task_array[i].f);
            direct.setFData(task_array[i].f_data);
            direct.setN(task_array[i].n);
            direct.setAB(task_array[i].A, task_array[i].B);
            direct.setMaxFeval(task_array[i].max_feval);
            direct.setMaxIter(task_array[i].max_iter);
            direct.setMagicEps(task_array[i].magic_eps);
            direct.setVolumeReltol(task_array[i].volume_reltol);
            direct.setSigmaReltol(task_array[i].sigma_reltol);
            direct.setLogfile(task_array[i].logfile);
            direct.setAlghorithm(task_array[i].algorithm);

            direct.solve(X, minf);

            cout << "Function: " << task_array[i].name << endl;
            cout << "Dimension = " << task_array[i].n << endl;
            cout << "[A; B] = [(" << task_array[i].A[0] << ", " << task_array[i].A[1] << "); (" << 
                                     task_array[i].B[0] << ", " << task_array[i].B[1] << ")]" << endl;
            cout << "X* = (" << task_array[i].X_opt[0] << ", " << task_array[i].X_opt[1] << ")" << endl;
            cout << "f(X*) = " << task_array[i].f(task_array[i].n, task_array[i].X_opt.data(),
                                                  &flag_tmp, &f_data_tmp) << endl;
            cout << "Parameters for method:" << endl;
            cout << "max_feval = " << max_feval << " max_iter = " << max_iter << endl;
            cout << "magic_eps = " << magic_eps << endl;
            cout << "volume_reltol = " << volume_reltol << " sigma_reltol = " << sigma_reltol << endl;
            cout << "Type of algorithm: " << ((task_array[i].algorithm == DIRECT_ORIGINAL) ? "DIRECT_ORIGINAL" : "DIRECT_GABLONSKY") << endl;
            cout << "Trials result:" << endl;
            cout << "Number of trials = " << count_trials << endl;
            cout << "X = (" << X[0] << ", " << X[1] << ")" << endl;
            cout << "f(X) = " << minf << endl;
            cout << "||X* - X|| = " << sqrt((task_array[i].X_opt[0] - X[0]) * (task_array[i].X_opt[0] - X[0]) + 
                                            (task_array[i].X_opt[1] - X[1]) * (task_array[i].X_opt[1] - X[1])) << endl;
            cout << "|f(X*) - f(X)| = " << abs(task_array[i].f(task_array[i].n, task_array[i].X_opt.data(), 
                                               &flag_tmp, &f_data_tmp) - minf) << endl;
            cout << endl;
        }
    }

#if defined(_MSC_VER)
    cin.get();
#endif

    return 0;
}
