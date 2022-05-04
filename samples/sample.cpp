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

using namespace std;

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
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f2(double x, int j) {
    switch (j) {
        case 1: return x * x - 0.05;
        case 2: return -x + 0.1;
        case 3: return 5.0 * x * x + 3.0 * x - 1.0;
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f3(double x, int j) {
    switch (j) {
        case 1: return sin(x);
        default: return numeric_limits<double>::quiet_NaN();
    }
}

int main() {
    ofstream ofstr("output_data_sample.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    vector<trial> trial_vec;

    double x_min, x_min_true;
    double a = 2.0, b = 8.0, eps = 0.001, r = 3.0;
    int m = 1, count;

    imgo_method imgo(f1, m, a, b, r, 0.0, eps);

    x_min_true = 2.0 * M_PI;

    x_min = imgo.solve(count);

    cout << "f1(x) = -2.0 * x + 3.0\n";
    cout << "g1(x) = sin(x)\n";
    print_result(f1, m, a, b, x_min_true, x_min, count);

    imgo.getTrialPoints(trial_vec);
    ofstr << a << " " << b << " " << m << " " << x_min << " " << x_min_true << endl;
    for (int i = 0; i < trial_vec.size(); i++) {
        ofstr << trial_vec[i].x << " " << trial_vec[i].z << endl;
    }
    ofstr << endl;

    m = 2;
    a = -2.0; b = 2.0;
    imgo.setFunc(f2);
    imgo.setAB(a, b);
    imgo.setM(m);
    x_min_true = 0.1;

    x_min = imgo.solve(count);

    cout << "f2(x) = 5.0 * x * x + 3.0 * x - 1.0\n";
    cout << "g1(x) = x * x - 0.05\n";
    cout << "g2(x) = -x + 0.1\n";
    print_result(f2, m, a, b, x_min_true, x_min, count);

    imgo.getTrialPoints(trial_vec);
    ofstr << a << " " << b << " " << m << " " << x_min << " " << x_min_true << endl;
    for (int i = 0; i < trial_vec.size(); i++) {
        ofstr << trial_vec[i].x << " " << trial_vec[i].z << endl;
    }
    ofstr << endl;

    m = 0;
    a = -4.0; b = 4.0;
    imgo.setFunc(f3);
    imgo.setAB(a, b);
    imgo.setM(m);
    x_min_true = -M_PI / 2.0;

    x_min = imgo.solve(count);

    cout << "f3(x) = sin(x)\n";
    print_result(f3, m, a, b, x_min_true, x_min, count);

    imgo.getTrialPoints(trial_vec);
    ofstr << a << " " << b << " " << m << " " << x_min << " " << x_min_true << endl;
    for (int i = 0; i < trial_vec.size(); i++) {
        ofstr << trial_vec[i].x << " " << trial_vec[i].z << endl;
    }
    ofstr << endl;

    ofstr.close();
#if defined( _MSC_VER )
    cin.get();
#endif
    return 0;
}
