#include<iostream>
#include<fstream>
#include<limits>
#if defined( _MSC_VER )
    #define _USE_MATH_DEFINES
    #include<math.h>
#else
    #include<cmath>
#endif

#include<imgo.h>
#include<task.h>

using std::ofstream;
using std::cerr;
using std::cout;
using std::endl;

void print_result(double (*f)(double, int), double m, double a, double b, double x_star, double x, int count) {
    cout << "[a; b] = [" << a << "; " << b << "]"<< endl;
    cout << "X* = " << x_star << endl;
    cout << "X = " << x << endl;
    cout << "|X* - X| = " << abs(x_star - x) << endl;
    cout << "|f(X*) - f(X)| = " << abs(f(x_star, m + 1) - f(x, m + 1)) << endl;
    cout << "Number of trials = " << count << endl;
    cout << endl;
}

double f1(double x, int j) {
    switch(j) {
        case 1: return sin(x);
        case 2: return -2.0 * x + 3.0;
        default: return std::numeric_limits<double>::quiet_NaN();
    }
}

double f2(double x, int j) {
    switch (j) {
        case 1: return x * x - 0.05;
        case 2: return -x + 0.1;
        case 3: return 5.0 * x * x + 3.0 * x - 1.0;
        default: return std::numeric_limits<double>::quiet_NaN();
    }
}

double f3(double x, int j) {
    switch (j) {
        case 1: return sin(x);
        default: return std::numeric_limits<double>::quiet_NaN();
    }
}

int main() {
    ofstream ofstr("output_data_sample.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    vector<trial_constr> trial_vec;

    double eps = 0.001, r = 3.0, d = 0.0, x;
    int count, Nmax = 1000;
    Stop stop = ACCURACY;

    vector<task_imgo> task_array = { task_imgo(f1, "f1(x)", 1, 2.0, 8.0, 2.0 * M_PI, eps, Nmax, r, d, stop),
                                     task_imgo(f2, "f2(x)", 2, -2.0, 2.0, 0.1, eps, Nmax, r, d, stop),
                                     task_imgo(f3, "f3(x)", 0, -4.0, 4.0, -M_PI / 2.0, eps ,Nmax, r, d, stop) };

    imgo_method imgo(nullptr);

    for (int i = 0; i < task_array.size(); i++) {
        imgo.setF(task_array[i].f);
        imgo.setM(task_array[i].m);
        imgo.setAB(task_array[i].A[0], task_array[i].B[0]);
        imgo.setEps(task_array[i].eps);
        imgo.setNmax(task_array[i].Nmax);
        imgo.setR(task_array[i].r);
        imgo.setD(task_array[i].d);

        imgo.solve(count, x, task_array[i].stop);

        cout << "Function: " << task_array[i].name << endl;
        cout << "[a; b] = [" << task_array[i].A[0] << "; " << task_array[i].B[0] << "]"<< endl;
        cout << "X* = " << task_array[i].X_opt[0] << endl;
        cout << "X = " << x << endl;
        cout << "|X* - X| = " << std::abs(task_array[i].X_opt[0] - x) << endl;
        cout << "|f(X*) - f(X)| = " << std::abs(task_array[i].f(task_array[i].X_opt[0], task_array[i].m + 1) - 
                                           task_array[i].f(x, task_array[i].m + 1)) << endl;
        cout << "Number of trials = " << count << endl;
        cout << endl;

        imgo.getTrialPoints(trial_vec);
        ofstr << task_array[i].A[0] << " " << task_array[i].B[0] << " " << task_array[i].m 
              << " " << x << " " << task_array[i].X_opt[0] << endl;
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
