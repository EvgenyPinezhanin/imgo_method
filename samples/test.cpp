#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#if defined( _MSC_VER )
    #define _USE_MATH_DEFINES
    #include<math.h>
#else
    #include<cmath>
#endif
#include<imgo.h>

using namespace std;

struct task {
    double (*f)(double, int);
    double a, b;
    double x_min;
    int m;

    task(double (*_f)(double, int), double _a, double _b, double _x_min, int _m)
        : f(_f), a(_a), b(_b), x_min(_x_min), m(_m) {}
};

double f1(double x, int j) {
    switch (j) {
        case 1: return exp(-sin(3.0 * x)) - 1.0 / 10.0 * pow(x - 1.0 / 2.0, 2.0) - 1.0;
        case 2: return -13.0 / 6.0 * x + sin(13.0 / 4.0 * (2.0 * x + 5.0)) - 53.0 / 12.0;
    }
    return -1.0;
}

double f2(double x, int j) {
    switch (j) {
        case 1: return 1.0 / 20.0 - exp(-2.0 / 5.0 * (x + 5.0)) * sin(4.0 / 5.0 * M_PI * (x + 5.0));
        case 2: return (11.0 * x * x - 10.0 * x + 21.0) / (2.0 * (x * x + 1));
    }
    return -1.0;
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
    }
    return -1.0;
}

// Error???
double f5(double x, int j) {
    switch (j) {
        case 1: return 17.0 / 25.0 - 2.0 / 29763.233 * (-1.0 / 6.0 * pow(x, 6) + 52.0 / 25.0 * pow(x, 5) - 39.0 / 80.0 * pow(x, 4) -
            71.0 / 10.0 * x * x * x + 79.0 / 20.0 * x * x + x - 1.0 / 10.0);
        case 2: return -14.0 / 125.0 * (3.0 * x - 8.0) * sin(252.0 / 125.0 * (x + 3.0 / 2.0)) - 1.0 / 2.0;
        case 3: return sin(0.423531 * x + 3.13531) + sin(10.0 / 3.0 * (0.423531 * x + 3.13531)) + log(0.423531 * x + 3.13531) + 0.36634 - 0.355766 * x;
    }
    return -1.0;
}

double f6(double x, int j) {
    switch (j) {
        case 1: return 40.0 * (cos(4.0 * x) * (x - sin(x)) * exp( -(x * x) / 2.0));
        case 2: return 2.0 / 25.0 * (x + 4.0) - sin(12.0 / 5.0 * (x + 4.0));
        case 3: return -7.0 / 40.0 * (3.0 * x + 4.0) * sin(63.0 / 20.0 * (x + 4.0));
    }
    return -1.0;
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
    }
    return -1.0;
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
    }
    return -1.0;
}

int main() {
    ofstream ofstr("output_data_test.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    vector<trial> trial_vec;

    double eps = 0.000001;
    double r = 2.0; // > 1

    double x_min;
    int n;

    vector<task> task_arr = {{f1, -2.5, 1.5, 1.05738, 1},
                             {f2, -5.0, 5.0, 1.016, 1},
                             {f4, 0.0, 4.0, 2.45956, 2},
                             {f5, -1.5, 11.0, 8.85725, 2},
                             {f6, -4.0, 4.0, 2.32396, 2},
                             {f8, -2.5, 1.5, -1.12724, 3},
                             {f9, 0.0, 14.0, 4.0, 3}};


    imgo_method imgo(&f1, 0, 0.0, 0.0, eps, r, 0.05);

    for (int i = 0; i < task_arr.size(); i++) {
        imgo.setFunc(task_arr[i].f);
        imgo.setA(task_arr[i].a);
        imgo.setB(task_arr[i].b);
        imgo.setM(task_arr[i].m);

        x_min = imgo.solve(n);

        cout << "f" << i + 1 << "(x)\n";
        cout << "[a; b] = [" << task_arr[i].a << "; " << task_arr[i].b << "]"<< endl;
        cout << "x_min = " << setprecision(12) << x_min << endl;
        cout << "Number of trials = " << n << endl;
        cout << "|Error rate| = " << abs(x_min - task_arr[i].x_min) << endl;
        cout << "|f" << i + 1 << "(x_min) - f" << i + 1 << "(x_opt)| = " <<
            abs(task_arr[i].f(x_min, task_arr[i].m + 1) - task_arr[i].f(task_arr[i].x_min, task_arr[i].m + 1)) << endl;
        cout << endl;

        imgo.getTrialPoints(trial_vec);
        ofstr << task_arr[i].a << " " << task_arr[i].b << " " << task_arr[i].m 
              << " " << x_min << " " << task_arr[i].x_min << endl;
        for (int j = 0; j < trial_vec.size(); j++) {
            ofstr << trial_vec[j].x << " " << trial_vec[j].z << endl;
        }
        ofstr << endl;
    }

    ofstr.close();
    cin.get();
    return 0;
}