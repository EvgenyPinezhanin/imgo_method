#ifndef TASK_H
#define TASK_H

#include <cstdio>
#include <string>
#include <vector>
#include <limits>

#include <opt_method.h>
#include <IGeneralOptProblem.hpp>
#include <IGeneralOptProblemFamily.hpp>
#include <IOptProblemFamily.hpp>
#include <IConstrainedOptProblemFamily.hpp>
#include <direct_method.h>

using namespace std;

struct task {
    string name;
    int n;
    vector<double> A, B, X_opt, L;
    double eps;
    int maxIters, maxEvals;

    bool used;

    task(string _name, int _n, vector<double> _A, vector<double> _B, vector<double> _X_opt, vector<double> _L, 
         double _eps, int _maxIters, int _maxEvals, bool _used = true) : name(_name), n(_n), A(_A), B(_B), X_opt(_X_opt), 
         L(_L), eps(_eps), maxIters(_maxIters), maxEvals(_maxEvals), used(_used) {};
};

struct task_gsa : public task {
    double (*f)(double);
    double r;

    task_gsa(double (*_f)(double), string _name, double _a, double _b, double _x_opt, double _L, double _eps, int _maxIters, 
             int _maxEvals, double _r, bool _used = true) : task(_name, 1, vector<double>{_a}, vector<double>{_b}, 
             vector<double>{_x_opt}, vector<double>{_L}, _eps, _maxIters, _maxEvals, _used), f(_f), r(_r) {};
};

struct task_imgo : public task {
    double (*f)(double, int);
    int numberConstraints;
    double r, d;

    task_imgo(double (*_f)(double, int), string _name, int _numberConstraints, double _a, double _b, double _x_opt,
              vector<double> _L, double _eps, int _maxIters, int _maxEvals, double _r, double _d, bool _used = true)
              : task(_name, 1, vector<double>{_a}, vector<double>{_b}, vector<double>{_x_opt}, _L, _eps, _maxIters,
              _maxEvals, _used), f(_f), numberConstraints(_numberConstraints), r(_r), d(_d) {};
};

struct task_mggsa : public task {
    double (*f)(vector<double>, int);
    int numberConstraints;
    double r, d;
    int den, key;

    task_mggsa(double (*_f)(vector<double>, int), string _name, int _n, int _numberConstraints, vector<double> _A,
               vector<double> _B, vector<double> _X_opt, vector<double> _L, double _eps, int _maxIters, int _maxEvals,
               double _r, double _d, int _den, int _key, bool _used = true) : task(_name, _n, _A, _B, _X_opt, _L, _eps,
               _maxIters, _maxEvals, _used), f(_f), numberConstraints(_numberConstraints), r(_r), d(_d), den(_den),
               key(_key) {};
};

enum class type_constraned { CONSTR, NONCONSTR };

struct problem {
    string name;
    type_constraned type;
    string short_name;

    bool used;

    problem(string _name, type_constraned _type, string _short_name = "", bool _used = true) 
            : name(_name), short_name(_short_name), type(_type), used(_used) {};
};

struct problem_single : public problem {
    IGeneralOptProblem *optProblem;

    problem_single(string _name, IGeneralOptProblem *_optProblem, type_constraned _type, string _short_name = "", bool _used = true) 
                   : problem(_name, _type, _short_name, used), optProblem(_optProblem) {};
};

struct problem_family : public problem {
    IGeneralOptProblemFamily *optProblemFamily;

    problem_family(string _name, IGeneralOptProblemFamily *_optProblemFamily, type_constraned _type, string _short_name = "", bool _used = true) 
                   : problem(_name, _type, _short_name, used), optProblemFamily(_optProblemFamily) {};
};

class Functor {
public:
    Functor() {};

    virtual double operator() (double x, int j) = 0;
    virtual double operator() (vector<double> x, int j) = 0; 
};

class functor_single : public Functor {
public:
    IOptProblem *opt_problem;

    functor_single(IOptProblem *_opt_problem = nullptr) 
                   : Functor(), opt_problem(_opt_problem) {};

    double operator() (double x, int j) override {
        switch (j) {
            case 1: return opt_problem->ComputeFunction(vector<double>{x});
            default: return numeric_limits<double>::quiet_NaN();
        }
    }

    double operator() (vector<double> x, int j) override {
        switch (j) {
            case 1: return opt_problem->ComputeFunction(x);
            default: return numeric_limits<double>::quiet_NaN();
        }
    }
};

class functor_single_constr : public Functor {
public:
    IConstrainedOptProblem *constr_opt_problem;

    functor_single_constr(IConstrainedOptProblem *_constr_opt_problem = nullptr) 
        : Functor(), constr_opt_problem(_constr_opt_problem) {}

    double operator() (double x, int j) override {
        int constr = constr_opt_problem->GetConstraintsNumber();
        if (j >= 1 && j <= constr) {
            return constr_opt_problem->ComputeConstraint(j - 1, vector<double>{x});
        } else if (j == constr + 1) {
            return constr_opt_problem->ComputeFunction(vector<double>{x});
        } else {
            return numeric_limits<double>::quiet_NaN();
        }
    }

    double operator() (vector<double> x, int j) override {
        int constr = constr_opt_problem->GetConstraintsNumber();
        if (j >= 1 && j <= constr) {
            return constr_opt_problem->ComputeConstraint(j - 1, x);
        } else if (j == constr + 1) {
            return constr_opt_problem->ComputeFunction(x);
        } else {
            return numeric_limits<double>::quiet_NaN();
        }
    }
};

class functor_family : public Functor {
public:
    IOptProblemFamily *opt_problem_family;
    int current_func;

    functor_family(IOptProblemFamily *_opt_problem_family = nullptr) 
        : Functor(), opt_problem_family(_opt_problem_family), current_func(0) {};

    double operator() (double x, int j) override {
        switch (j) {
            case 1: return (*opt_problem_family)[current_func]->ComputeFunction(vector<double>{x});
            default: return numeric_limits<double>::quiet_NaN();
        }
    }

    double operator() (vector<double> x, int j) override {
        switch (j) {
            case 1: return (*opt_problem_family)[current_func]->ComputeFunction(x);
            default: return numeric_limits<double>::quiet_NaN();
        }
    }
};

class functor_family_constr : public Functor {
public:
    IConstrainedOptProblemFamily *constr_opt_problem_family;
    int current_func;

    functor_family_constr(IConstrainedOptProblemFamily *_constr_opt_problem_family = nullptr) 
        : Functor(), constr_opt_problem_family(_constr_opt_problem_family), current_func(0) {}

    double operator() (double x, int j) {
        int constr = (*constr_opt_problem_family)[current_func]->GetConstraintsNumber();
        if (j >= 1 && j <= constr) {
            return (*constr_opt_problem_family)[current_func]->ComputeConstraint(j - 1, vector<double>{x});
        } else if (j == constr + 1) {
            return (*constr_opt_problem_family)[current_func]->ComputeFunction(vector<double>{x});
        } else {
            return numeric_limits<double>::quiet_NaN();
        }
    }

    double operator() (vector<double> x, int j) {
        int constr = (*constr_opt_problem_family)[current_func]->GetConstraintsNumber();
        if (j >= 1 && j <= constr) {
            return (*constr_opt_problem_family)[current_func]->ComputeConstraint(j - 1, x);
        } else if (j == constr + 1) {
            return (*constr_opt_problem_family)[current_func]->ComputeFunction(x);
        } else {
            return numeric_limits<double>::quiet_NaN();
        }
    }
};

struct data_direct {
    int count_evals;
    vector<vector<double>> points;
};

struct data_direct_oper_character : public data_direct {
    Functor *functor;
    type_constraned type;
    bool converge;
    int min_count_evals;
    double eps;
};

struct task_direct {
    string name;

    direct_objective_func f;
    void *f_data;
    int n;
    vector<double> A, B, X_opt;

    int max_feval, max_iter;
    double magic_eps, magic_eps_abs;
    double volume_reltol, sigma_reltol;

    FILE *logfile;
    direct_algorithm algorithm;

    bool used;

    task_direct(string _name, direct_objective_func _f, void *_f_data, int _n, vector<double> _A, vector<double> _B, 
                vector<double> _X_opt, int _max_feval, int _max_iter, double _magic_eps, double _magic_eps_abs,
                double _volume_reltol, double _sigma_reltol, FILE *_logfile, direct_algorithm _algorithm, bool _used = true) 
                : name(_name), f(_f), f_data(_f_data), n(_n), A(_A), B(_B), X_opt(_X_opt), max_feval(_max_feval), 
                max_iter(_max_iter), magic_eps(_magic_eps), magic_eps_abs(_magic_eps_abs), volume_reltol(_volume_reltol), 
                sigma_reltol(_sigma_reltol), logfile(_logfile), algorithm(_algorithm), used(_used) {};
};

#endif // TASK_H
