#ifndef IMGO_H
#define IMGO_H

#include<vector>
#include<algorithm>
#include<cmath>

using namespace std;

struct trial {
    double x, z;
    size_t nu;
    //int H_number;

    trial(double _x = 0.0, double _z = 0.0, int _nu = 0) //, int _H_number = -1) 
        : x(_x), z(_z), nu(_nu) {} //, H_number(_H_number) {};
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

    //test
    vector<trial> trial_points_test;
    //test

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

    //test
    void getTrialPointsTest(vector<trial>& trial_vec);
    //test
    
    double solve(int &n);
    bool solve_test(double x_opt, int k);
};

/*
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
*/
#endif