#ifndef IMGO_H
#define IMGO_H

#include<vector>
#include<algorithm>
#include<cmath>

using namespace std;

struct trial {
    double x, z;
    int v;
    int H_number;

    trial(double _x = 0.0, double _z = 0.0, int _v = 0, int _H_number = -1) 
        : x(_x), z(_z), v(_v), H_number(_H_number) {};
};

class imgo_method {
protected:
    double a, b;
    double eps;
    double r;
    int m;
    double (*f)(double, int);

    double d;

    vector<trial> a_inf;
    vector<vector<trial>> I;
    vector<double> z_v;
    vector<double> mu;

    void addInSort(vector<trial> &vec, trial tr);
    double searchMinX();
    trial trial_func(double x);
    double selectNewPoint(int &t);

public:
    imgo_method(double (*_f)(double, int), int _m, double _a, double _b, double _eps, double _r, double _d = 0.01);
    
    void setFunc(double (*_f)(double, int));
    void setA(double _a);
    void setB(double _b);
    void setM(int _m);
    
    double solve(int &n);
};

class imgo_method_adaptive : public imgo_method {
private:
    vector<vector<int>> H;
    vector<int> H_default;

    void changeH(vector<int> &H, int p);
    void changeH(vector<int> &H, int p, int q);
    int getNumberOfIndex(double v, trial tr);
    void trial_func(trial &tr);
    trial selectNewPoint(int &t);

public:
    imgo_method_adaptive(double (*_f)(double, int), int _m, double _a, double _b, double _eps, double _r, double _d = 0.01);

    void setM(int _m);
    
    double solve(int &n);
};

#endif