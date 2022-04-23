#ifndef TASK_H
#define TASK_H

#include <string>
#include <vector>
#include <imgo.h>

using namespace std;

struct ITask {
    string name;
    int n, m;
    vector<double> A, B;
    vector<double> X_opt;
    double eps, r, d;
    int Nmax;
    Stop stop;
    bool used;

    ITask(string _name, int _n, int _m, vector<double> _A, vector<double> _B, vector<double> _X_opt, 
    double _eps, int _Nmax, double _r, double _d, Stop _stop, bool _used = true)
        : name(_name), n(_n), m(_m), A(_A), B(_B), X_opt(_X_opt), eps(_eps), Nmax(_Nmax), r(_r), d(_d), stop(_stop), used(_used) {};
};

struct Task : public ITask {
    double (*f)(double, int);

    Task(double (*_f)(double, int), string _name, int _m, double _a, double _b, double _x, 
    double _eps, int _Nmax, double _r, double _d, Stop _stop, bool _used = true)
        : ITask(_name, 1, _m, vector<double>{_a}, vector<double>{_b}, vector<double>{_x}, _eps, _Nmax, _r, _d, _stop, _used), f(_f) {};
};

struct Task_peano : public ITask {
    double (*f)(vector<double>, int);
    int den, key;

    Task_peano(double (*_f)(vector<double>, int), string _name, int _n, int _m, vector<double> _A, vector<double> _B, 
    vector<double> _X_opt, double _eps, int _Nmax, double _r, double _d, int _den, int _key, Stop _stop, bool _used = true)
        : ITask(_name, _n, _m, _A, _B, _X_opt, _eps, _Nmax, _r, _d, _stop, _used), f(_f), den(_den), key(_key) {};
};

#endif // TASK_H
