#ifndef IMGO_H
#define IMGO_H

#include<vector>
#include<algorithm>
#include<cmath>

using namespace std;

struct trial {
    double x, z;
    size_t nu;

    trial(double _x = 0.0, double _z = 0.0, int _nu = 0) : x(_x), z(_z), nu(_nu) {}
};

class imgo_method {
protected:
    double (*f)(double, int);
    size_t m;
    double a, b;
    double eps;
    double r;
    double d;

    vector<trial> trial_points;
    vector<vector<trial>> I;
    vector<double> z_star;
    vector<double> mu;

    void addInSort(vector<trial> &vec, trial tr);
    double searchMinX();
    trial newTrial(double x);
    double selectNewPoint(size_t &t);

public:
    imgo_method(double (*_f)(double, int), size_t _m = 0, double _a = 0.0, double _b = 10.0, double _eps = 0.0001, double _r = 2.0, double _d = 0.01);
    
    void setFunc(double (*_f)(double, int));
    void setA(double _a);
    void setB(double _b);
    void setM(size_t _m);
    void setEps(double _eps);

    void getTrialPoints(vector<trial> &trial_vec);
    
    double solve(int &n);
    bool solve_test(double x_opt, int k);
};

#endif