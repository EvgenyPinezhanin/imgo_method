#if defined( _MSC_VER )
    #define _CRT_SECURE_NO_WARNINGS    
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include <direct_method.h>
#include <task.h>

using namespace std;

const int test_func_number = 0; // 0 - f1, 1 - f2

double f1(int n, const double *X, int *undefined_flag, void *data) {
    data_direct *f_data = static_cast<data_direct*>(data);
    f_data->count_evals++;
    f_data->points.push_back(vector<double>{X[0], X[1]});

    bool is_constr = true;
    vector<double> constr(3);
    constr[0] = 0.01 * (pow(X[0] - 2.2, 2) + pow(X[1] - 1.2, 2) - 2.25);
    constr[1] = 100.0 * (1.0 - pow(X[0] - 2.0, 2) / 1.44 - pow(0.5 * X[1], 2));
    constr[2] = 10.0 * (X[1] - 1.5 - 1.5 * sin(6.283 * (X[0] - 1.75)));
    for (int i = 0; i < 3; i++) {
        if (constr[i] > 0.0) is_constr = false;
    }
    if (!is_constr) *undefined_flag = 1;

    return -1.5 * X[0] * X[0] * exp(1.0 - X[0] * X[0] - 20.25 * pow((X[0] - X[1]), 2)) -
            pow(0.5 * (X[0] - 1.0) * (X[1] - 1.0), 4) * exp(2.0 - pow(0.5 * (X[0] - 1.0), 4) -
            pow(X[1] - 1.0, 4));
}

const double C[20] = {75.1963666677,-3.8112755343,0.1269366345,-0.0020567665,0.000010345,
                      -6.8306567631,0.0302344793,-0.0012813448,0.0000352559,-0.0000002266,
                      0.2564581253,-0.0034604030,0.0000135139,-28.1064434908,-0.0000052375,
                      -0.0000000063,0.0000000007,0.0003405462,-0.0000016638,-2.8673112392 };

double f2(int n, const double *X, int *undefined_flag, void *data) {
    data_direct *f_data = static_cast<data_direct*>(data);
    f_data->count_evals++;
    f_data->points.push_back(vector<double>{X[0], X[1]});

    bool is_constr = true;
    vector<double> constr(4);
    constr[0] = 450.0 - X[0] * X[1];
    constr[1] = (0.1 * X[0] - 1.0) * (0.1 * X[0] - 1.0) - X[1];
    constr[2] = 8.0 * (X[0] - 40.0) - (X[1] - 30.0) * (X[1] - 55.0);
    constr[3] = X[1] + (X[0] - 35.0) * (X[0] - 30.0) / 125.0 - 80.0;
    for (int i = 0; i < 4; i++) {
        if (constr[i] > 0.0) is_constr = false;
    }
    if (!is_constr) *undefined_flag = 1;

    return -(C[0] + C[1] * X[0] + C[2] * X[0] * X[0] + C[3] * pow(X[0], 3) + C[4] * pow(X[0], 4) + C[5] * X[1] +
             C[6] * X[0] * X[1] + C[7] * X[0] * X[0] * X[1] + C[8] * pow(X[0], 3) * X[1] + C[9] * pow(X[0], 4) * X[1] +
             C[10] * X[1] * X[1] + C[11] * pow(X[1], 3) + C[12] * pow(X[1], 4) + C[13] / (X[1] + 1) + C[14] * X[0] * X[0] * X[1] *X[1] +
             C[15] * pow(X[0], 3) * X[1] * X[1] + C[16] * pow(X[0], 3) * pow(X[1], 3) + C[17] * X[0] * X[1] * X[1] +
             C[18] * X[0] * pow(X[1], 3) + C[19] * exp(0.0005 * X[0] * X[1]));
}

int main() {
    ofstream ofstr("output_data/direct_test.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstr_opt("output_data/direct_test_opt.txt");
    if (!ofstr_opt.is_open()) cerr << "File opening error\n";

    data_direct f_data;
    int n = 2;
    int max_feval = 1000, max_iter = 1000;
    double magic_eps = 1.0e-4;
    double method_accuracy = 1.0e-5;
    double volume_reltol = method_accuracy / sqrt(n);
    double sigma_reltol = 0.0;
    direct_algorithm algorithm = DIRECT_ORIGINAL;
    vector<double> X;
    double minf;

    vector<task_direct> task_array = { task_direct("f1", f1, &f_data, n, vector<double>{0.0, -1.0}, vector<double>{4.0, 3.0},
                                                   vector<double>{0.942, 0.944}, max_feval, max_iter, magic_eps, volume_reltol,
                                                   sigma_reltol, nullptr, algorithm),
                                       task_direct("f2", f2, &f_data, n, vector<double>{0.0, 0.0}, vector<double>{80.0, 80.0},
                                                   vector<double>{77.489, 63.858}, max_feval, max_iter, magic_eps, volume_reltol,
                                                   sigma_reltol, nullptr, algorithm) };

    direct_method direct;

    int flag_tmp;
    data_direct *data, data_tmp;
    for (int i = 0; i < task_array.size(); i++) {
        if (task_array[i].used) {
            data = static_cast<data_direct*>(task_array[i].f_data);
            data->count_evals = 0;
            data->points.clear();

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
                                                  &flag_tmp, &data_tmp) << endl;
            cout << "Parameters for method:" << endl;
            cout << "max_feval = " << task_array[i].max_feval << " max_iter = " << task_array[i].max_iter << endl;
            cout << "magic_eps = " << task_array[i].magic_eps << endl;
            cout << "volume_reltol = " << task_array[i].volume_reltol << " sigma_reltol = " << task_array[i].sigma_reltol << endl;
            cout << "Type of algorithm: " << ((task_array[i].algorithm == DIRECT_ORIGINAL) ? "DIRECT_ORIGINAL" : "DIRECT_GABLONSKY") << endl;
            cout << "Trials result:" << endl;
            cout << "Number of trials = " << data->count_evals << endl;
            cout << "X = (" << X[0] << ", " << X[1] << ")" << endl;
            cout << "f(X) = " << minf << endl;
            cout << "||X* - X|| = " << sqrt((task_array[i].X_opt[0] - X[0]) * (task_array[i].X_opt[0] - X[0]) + 
                                            (task_array[i].X_opt[1] - X[1]) * (task_array[i].X_opt[1] - X[1])) << endl;
            cout << "|f(X*) - f(X)| = " << abs(task_array[i].f(task_array[i].n, task_array[i].X_opt.data(), 
                                               &flag_tmp, &data_tmp) - minf) << endl;
            cout << endl;

            // Saving points for plotting
            ofstr << task_array[i].X_opt[0] << " " << task_array[i].X_opt[1] << " " << 
                     task_array[i].f(task_array[i].n, task_array[i].X_opt.data(), &flag_tmp, &data_tmp) << endl;
            ofstr << endl << endl;
            ofstr << X[0] << " " << X[1] << " " << task_array[i].f(task_array[i].n, X.data(), &flag_tmp, &data_tmp) << endl;
            ofstr << endl << endl;
            for (int j = 0; j < data->points.size(); j++) {
                ofstr << data->points[j][0] << " " << data->points[j][1] << " " << 
                         task_array[i].f(task_array[i].n, data->points[j].data(), &flag_tmp, &data_tmp) << endl;
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
    error = system("chmod +x scripts/direct_test.gp");
    if (error != 0) {
        cerr << "Error chmod" << endl;
    }
#endif

    char str[100];
    sprintf(str, "gnuplot -c scripts/direct_test.gp %d", test_func_number);
    error = system(str);
    if (error != 0) {
        cerr << "Error gnuplot" << endl;
    }

#if defined(_MSC_VER)
    cin.get();
#endif

    return 0;
}
