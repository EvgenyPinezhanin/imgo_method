#ifndef TASK_H
#define TASK_H

#include <string>
#include <vector>
#include <limits>

#include <opt_method.h>
#include <IGeneralOptProblem.hpp>
#include <IGeneralOptProblemFamily.hpp>
#include <IOptProblemFamily.hpp>
#include <IConstrainedOptProblemFamily.hpp>

using namespace std;

struct task {
    string name;
    int n;
    vector<double> A, B, X_opt, L;
    double eps;
    int Nmax;
    Stop stop;

    bool used;

    task(string _name, int _n, vector<double> _A, vector<double> _B, vector<double> _X_opt, vector<double> _L, 
         double _eps, int _Nmax, Stop _stop, bool _used = true) : name(_name), n(_n), A(_A), B(_B), X_opt(_X_opt), 
         L(_L), eps(_eps), Nmax(_Nmax), stop(_stop), used(_used) {};
};

struct task_gsa : public task {
    double (*f)(double);
    double r;

    task_gsa(double (*_f)(double), string _name, double _a, double _b, double _x_opt, double _L, double _eps, int _Nmax, 
             double _r, Stop _stop, bool _used = true) : task(_name, 1, vector<double>{_a}, vector<double>{_b}, 
             vector<double>{_x_opt}, vector<double>{_L}, _eps, _Nmax, _stop, _used), f(_f), r(_r) {};
};

struct task_imgo : public task {
    double (*f)(double, int);
    int m;
    double r, d;

    task_imgo(double (*_f)(double, int), string _name, int _m, double _a, double _b, double _x_opt, vector<double> _L,
              double _eps, int _Nmax, double _r, double _d, Stop _stop, bool _used = true) : task(_name, 1, vector<double>{_a}, 
              vector<double>{_b}, vector<double>{_x_opt}, _L, _eps, _Nmax, _stop, _used), f(_f), m(_m), r(_r), d(_d) {};
};

struct task_mggsa : public task {
    double (*f)(vector<double>, int);
    int m;
    double r, d;
    int den, key;

    task_mggsa(double (*_f)(vector<double>, int), string _name, int _n, int _m, vector<double> _A, vector<double> _B, 
               vector<double> _X_opt, vector<double> _L, double _eps, int _Nmax, double _r, double _d, int _den, int _key, 
               Stop _stop, bool _used = true) : task(_name, _n, _A, _B, _X_opt, _L, _eps, _Nmax, _stop, _used), f(_f), 
               m(_m), r(_r), d(_d), den(_den), key(_key) {};
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

struct functor_non_constr {
    functor_non_constr(IOptProblemFamily *_opt_problem_family = nullptr) 
        : opt_problem_family(_opt_problem_family), current_func(0) {};

    IOptProblemFamily *opt_problem_family;
    int current_func;

    double operator() (double x, int j) {
        switch (j) {
            case 1: return (*opt_problem_family)[current_func]->ComputeFunction(vector<double>{x});
            default: return numeric_limits<double>::quiet_NaN();
        }
    }

    double operator() (vector<double> x, int j) {
        switch (j) {
            case 1: return (*opt_problem_family)[current_func]->ComputeFunction(x);
            default: return numeric_limits<double>::quiet_NaN();
        }
    }
};

struct functor_constr {
    functor_constr(IConstrainedOptProblemFamily *_constr_opt_problem_family = nullptr) 
        : constr_opt_problem_family(_constr_opt_problem_family), current_func(0) {}

    IConstrainedOptProblemFamily *constr_opt_problem_family;
    int current_func;

    double operator() (double x, int j) {
        int constr = (*constr_opt_problem_family)[current_func]->GetConstraintsNumber();
        if (j >= 1 && j <= constr) {
            return (*constr_opt_problem_family)[current_func]->ComputeConstraint(j - 1, vector<double>(x));
        } else if (j - 1 == constr) {
            return (*constr_opt_problem_family)[current_func]->ComputeFunction(vector<double>(x));
        } else {
            return numeric_limits<double>::quiet_NaN();
        }
    }

    double operator() (vector<double> x, int j) {
        int constr = (*constr_opt_problem_family)[current_func]->GetConstraintsNumber();
        if (j >= 1 && j <= constr) {
            return (*constr_opt_problem_family)[current_func]->ComputeConstraint(j - 1, x);
        } else if (j - 1 == constr) {
            return (*constr_opt_problem_family)[current_func]->ComputeFunction(x);
        } else {
            return numeric_limits<double>::quiet_NaN();
        }
    }
};

#endif // TASK_H
