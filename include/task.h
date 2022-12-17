#ifndef TASK_H
#define TASK_H

#include <string>
#include <vector>

#include <opt_method.h>
#include <IGeneralOptProblem.hpp>
#include <IGeneralOptProblemFamily.hpp>

using namespace std;

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

    task_imgo(double (*_f)(double, int), string _name, int _m, double _a, double _b, double _x_opt, 
              double _eps, int _Nmax, double _r, double _d, Stop _stop, bool _used = true)
             : task(_name, 1, vector<double>{_a}, vector<double>{_b}, vector<double>{_x_opt}, 
                    _eps, _Nmax, _stop, _used), f(_f), m(_m), r(_r), d(_d) {};
};

struct task_mggsa : public task {
    double (*f)(vector<double>, int);
    int m;
    double r, d;
    int den, key;

    task_mggsa(double (*_f)(vector<double>, int), string _name, int _n, int _m, vector<double> _A, vector<double> _B, 
               vector<double> _X_opt, double _eps, int _Nmax, double _r, double _d, int _den, int _key, Stop _stop, bool _used = true)
              : task(_name, _n, _A, _B, _X_opt, _eps, _Nmax, _stop, _used), f(_f), m(_m), r(_r), d(_d), den(_den), key(_key) {};
};

enum class type_constraned { CONSTR, NONCONSTR };

struct class_problems {
    string name;
    string short_name;
    type_constraned type;

    bool used;

    class_problems(string _name, type_constraned _type, string _short_name = "def", bool _used = true) 
        : name(_name), short_name(_short_name), type(_type), used(_used) {};
};

struct class_problems_o : public class_problems {
    IGeneralOptProblem *problem;

    class_problems_o(string _name, IGeneralOptProblem *_problem, type_constraned _type, string _short_name = "def", bool _used = true) 
        : class_problems(_name, _type, _short_name, used), problem(_problem) {};
};

struct class_problems_os : public class_problems_o {
    double (*f)(double, int);

    class_problems_os(string _name, IGeneralOptProblem *_problem, type_constraned _type, double (*_f)(double, int), string _short_name = "def",
        bool _used = true) : class_problems_o(_name, _problem, _type, _short_name, _used), f(_f) {};
};

struct class_problems_om : public class_problems_o {
    double (*f)(vector<double>, int);

    class_problems_om(string _name, IGeneralOptProblem *_problem, type_constraned _type, double (*_f)(vector<double>, int), 
    string _short_name = "def", bool _used = true) : class_problems_o(_name, _problem, _type, _short_name, _used), f(_f) {};
};

struct class_problems_f : public class_problems {
    IGeneralOptProblemFamily *problem;

    class_problems_f(string _name, IGeneralOptProblemFamily *_problem, type_constraned _type, string _short_name = "def", bool _used = true) 
        : class_problems(_name, _type, _short_name, used), problem(_problem) {};
};

struct class_problems_fs : public class_problems_f {
    double (*f)(double, int);

    class_problems_fs(string _name, IGeneralOptProblemFamily *_problem, type_constraned _type, double (*_f)(double, int), 
        string _short_name = "def", bool _used = true) : class_problems_f(_name, _problem, _type, _short_name, _used), f(_f) {};
};

struct class_problems_fm : public class_problems_f {
    double (*f)(vector<double>, int);

    class_problems_fm(string _name, IGeneralOptProblemFamily *_problem, type_constraned _type, double (*_f)(vector<double>, int), 
    string _short_name = "def", bool _used = true) : class_problems_f(_name, _problem, _type, _short_name, _used), f(_f) {};
};

#endif // TASK_H
