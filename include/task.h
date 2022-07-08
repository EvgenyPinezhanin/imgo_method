#ifndef TASK_H
#define TASK_H

#include <string>
#include <vector>

#include <opt_method.h>

using std::string;
using std::vector;

struct task {
    string name;
    int n;
    vector<double> A, B;
    vector<double> X_opt;
    double eps;
    int Nmax;
    Stop stop;

    bool used;

    task(string _name, int _n, vector<double> _A, vector<double> _B, vector<double> _X_opt, double _eps, int _Nmax, Stop _stop, bool _used = true)
        : name(_name), n(_n), A(_A), B(_B), X_opt(_X_opt), eps(_eps), Nmax(_Nmax), stop(_stop), used(_used) {};
};

struct task_gsa : public task {
    double (*f)(double);
    double r;

    task_gsa(double (*_f)(double), string _name, double _a, double _b, double _x_opt, 
             double _eps, int _Nmax, double _r, Stop _stop, bool _used = true)
            : task(_name, 1, vector<double>{_a}, vector<double>{_b}, vector<double>{_x_opt}, _eps, _Nmax, _stop, _used), f(_f), r(_r) {};
};

struct task_imgo : public task {
    double (*f)(double, int);
    int m;
    double r, d;

    task_imgo(double (*_f)(double, int), string _name, int _m, vector<double> _A, vector<double> _B, 
              vector<double> _X_opt, double _eps, int _Nmax, double _r, double _d, Stop _stop, bool _used = true)
             : task(_name, 1, _A, _B, _X_opt, _eps, _Nmax, _stop, _used), f(_f), m(_m), r(_r), d(_d) {};
};

struct task_mggsa : public task {
    double (*f)(vector<double>, int);
    int m;
    double r, d;
    int den, key;

    task_mggsa(double (*_f)(vector<double>, int), string _name, int _n, int _m, vector<double> _A, vector<double> _B, 
               vector<double> _X_opt, double _eps, int _Nmax, double _r, double _d, int _den, int _key, Stop _stop, bool _used = true)
              : task(_name, _n, _A, _B, _X_opt, _eps, _Nmax, _stop, _used), f(_f), den(_den), key(_key) {};
};

#endif // TASK_H
