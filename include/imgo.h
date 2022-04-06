#ifndef IMGO_H
#define IMGO_H

#include <vector>
#include <algorithm>
#include <cmath>
#include <peano.h>

using namespace std;

enum Stop {ACCURACY, NUMBER, ACCURNUMBER}; // сделать ACCURNUMBER

struct trial {
    double x, z;
    size_t nu;

    trial(double _x = 0.0, double _z = 0.0, int _nu = 0) : x(_x), z(_z), nu(_nu) {}
};

class imgo_method {
protected:
    double (*f)(double, int);
    double (*f_md)(vector<double>, int);
    size_t m;
    int n, den, key;

    double a, b;
    vector<double> A, B;
    double eps, r, d;
    int Nmax;

    int M;

    vector<trial> trial_points;
    vector<vector<trial>> I;
    vector<double> z_star;
    vector<double> mu;

    void addInSort(vector<trial> &vec, trial tr);
    double searchMinX();
    trial newTrial(double x);
    trial newTrial_md(double x);
    double deltaX(double x1, double x2);
    double newPoint(size_t t);
    double selectNewPoint(size_t &t);

public:
    imgo_method(double (*_f)(double, int), size_t _m = 0, double _a = 0.0, double _b = 10.0, double _eps = 0.0001, double _r = 2.0, double _d = 0.01);
    imgo_method(double (*_f)(vector<double>, int), int _n, size_t _m = 0, vector<double> _A = vector<double>(), 
                vector<double> _B = vector<double>(), double _eps = 0.0001, double _r = 2.0, double _d = 0.01, int _den = 10, int _key = 1);
    
    void setFunc(double (*_f)(double, int));
    void setFunc(double (*_f_md)(vector<double>, int));
    void setA(double _a);
    void setA(vector<double> _A);
    void setB(double _b);
    void setB(vector<double> _B);
    void setM(size_t _m);
    void setEps(double _eps);
    void setR(double _r);
    void setD(double _d);
    void setN(int _n);
    void setDen(int _den);
    void setKey(int _key);


    void getTrialPoints(vector<trial> &trial_vec);
    void getPoints(vector<vector<double>> &points_vec);

    void y(double x, vector<double> &X);
    
    double solve(int &n, Stop stop = ACCURACY);
    void solve(int &n, vector<double> &X, Stop stop = ACCURACY);
    bool solve_test(double x_opt, int k);
    bool solve_test(vector<double> x_opt, int k);
};

#endif // IMGO_H
