#ifndef IMGO_H
#define IMGO_H

#include<vector>
#include<algorithm>
#include<cmath>

using namespace std;

const double eps = 1;

struct trial {
    double x, z;
    int v;

    trial(double _x, double _z, int _v) : x(_x), z(_z), v(_v) {};
};

class imgo_method {
private:
    double a, b;
    double eps;
    double r;
    int m;
    double (*f)(double, int);

    vector<trial> a_inf;
    vector<vector<trial>> I;
    vector<double> z_v;
    vector<double> mu;

    void addInSort(vector<trial> &vec, trial tr);
    double searchMinX();
    trial trial_func(double x);
    double selectNewPoint(int &t);

public:
    imgo_method(double (*_f)(double, int), int _m, double _a, double _b, double _eps, double _r);
    
    void setFunc(double (*_f)(double, int));
    void setA(double _a);
    void setB(double _b);
    void setM(int _m);
    
    double solve(int &n);
};

#endif