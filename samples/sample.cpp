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

    double a = 2.0, b = 8.0;
    double eps = 0.001;
    double r = 3;
    int m;

    double x_min;
    double x_min_true;
    int n;

    m = 1;
    imgo_method imgo(f1, m, a, b, r, 0.0,  eps);
    x_min = imgo.solve(n);
    x_min_true = 2.0 * M_PI;

    cout << "f1(x) = -2.0 * x + 3.0\n";
    cout << "g1(x) = sin(x)\n";
    cout << "[a; b] = [" << a << "; " << b << "]"<< endl;
    cout << "X_min_true = " << x_min_true << endl;
    cout << "X_min = " << x_min << endl;
    cout << "Number of trials = " << n << endl;
    cout << "|Error rate| = " << abs(x_min - x_min_true) << endl;
    cout << endl;

    imgo.getTrialPoints(trial_vec);
    ofstr << a << " " << b << " " << m << " " << x_min << " " << x_min_true << endl;
    for (int i = 0; i < trial_vec.size(); i++) {
        ofstr << trial_vec[i].x << " " << trial_vec[i].z << endl;
    }
    ofstr << endl;

    m = 2;
    a = -2.0; b = 2.0;
    imgo.setFunc(&f2);
    imgo.setA(a);
    imgo.setB(b);
    imgo.setM(m);
    x_min = imgo.solve(n);
    x_min_true = 0.1;

    cout << "f2(x) = 5.0 * x * x + 3.0 * x - 1.0\n";
    cout << "g1(x) = x * x - 0.05\n";
    cout << "g2(x) = -x + 0.1\n";
    cout << "[a; b] = [" << a << "; " << b << "]"<< endl;
    cout << "X_min = " << x_min << endl;
    cout << "Number of trials = " << n << endl;
    cout << "|Error rate| = " << abs(x_min - x_min_true) << endl;
    cout << endl;

    imgo.getTrialPoints(trial_vec);
    ofstr << a << " " << b << " " << m << " " << x_min << " " << x_min_true << endl;
    for (int i = 0; i < trial_vec.size(); i++) {
        ofstr << trial_vec[i].x << " " << trial_vec[i].z << endl;
    }
    ofstr << endl;

    m = 0;
    a = -4.0; b = 4.0;
    imgo.setFunc(&f3);
    imgo.setA(a);
    imgo.setB(b);
    imgo.setM(m);
    x_min = imgo.solve(n);
    x_min_true = -M_PI / 2.0;

    cout << "f3(x) = sin(x)\n";
    cout << "[a; b] = [" << a << "; " << b << "]" << endl;
    cout << "X_min = " << x_min << endl;
    cout << "Number of trials = " << n << endl;
    cout << "|Error rate| = " << abs(x_min - x_min_true) << endl;
    cout << endl;

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