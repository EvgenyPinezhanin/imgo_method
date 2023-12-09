#include <vector>
#include <string>

// #include <tasks/Task.h>
// #include <opt_problems/OneDimensionalConstrainedProblem.h>

using std::vector;
using std::string;


#ifndef STRONGIN_TASK_H_
#define STRONGIN_TASK_H_

// #include <base_structures/tasks/GeneralNumericalTask.h>

/* template <typename OptimizationProblemType>
struct StronginTask : public GeneralNumericalTask<OptimizationProblemType> {
    double r, constantEstimation;

    StronginTask(string _name, const OptimizationProblemType &_optTask, double _r, double _constantEstimation,
                 double _accuracy, int _maxTrials, int _maxFevals, bool _use)
                 : GeneralNumericalTask<OptimizationProblemType>(_name, _optTask, _accuracy, _maxTrials, _maxFevals, _use),
                 r(_r), constantEstimation(_constantEstimation) {};
}; */

#endif // STRONGIN_TASK_H_

/*

struct TaskMggsa : public Task {
    double (*f)(vector<double>, int);
    int numberConstraints;
    double r, d;
    int den, key, incr;

    TaskMggsa(double (*_f)(vector<double>, int), string _name, int _n, int _numberConstraints, vector<double> _A,
              vector<double> _B, vector<double> _XOpt, vector<double> _L, double _eps, int _maxTrials, int _maxFevals,
              double _r, double _d, int _den, int _key, int _incr, bool _used = true) : Task(_name, _n, _A, _B, _XOpt, _L, _eps,
              _maxTrials, _maxFevals, _used), f(_f), numberConstraints(_numberConstraints), r(_r), d(_d), den(_den),
              key(_key), incr(_incr) {};
};

enum class TypeConstraints { Constraints, NoConstraints };

struct Problem {
    string name;
    TypeConstraints type;
    string shortName;

    bool used;

    Problem(string _name, TypeConstraints _type, string _shortName = "", bool _used = true) 
            : name(_name), shortName(_shortName), type(_type), used(_used) {};
};

struct ProblemSingle : public Problem {
    IGeneralOptProblem *optProblem;

    ProblemSingle(string _name, IGeneralOptProblem *_optProblem, TypeConstraints _type, string _shortName = "", bool _used = true) 
                  : Problem(_name, _type, _shortName, used), optProblem(_optProblem) {};
};

class Functor {
public:
    Functor() {};

    virtual double operator() (double x, int j) = 0;
    virtual double operator() (vector<double> x, int j) = 0; 
};

class FunctorSingle : public Functor {
public:
    IOptProblem *optProblem;

    FunctorSingle(IOptProblem *_optProblem = nullptr) : Functor(), optProblem(_optProblem) {};

    double operator() (double x, int j) override {
        switch (j) {
            case 1: return optProblem->ComputeFunction(vector<double>{ x });
            default: return numeric_limits<double>::quiet_NaN();
        }
    }

    double operator() (vector<double> x, int j) override {
        switch (j) {
            case 1: return optProblem->ComputeFunction(x);
            default: return numeric_limits<double>::quiet_NaN();
        }
    }
};

class FunctorSingleConstrained : public Functor {
public:
    IConstrainedOptProblem *constrainedOptProblem;

    FunctorSingleConstrained(IConstrainedOptProblem *_constrainedOptProblem = nullptr) 
                             : Functor(), constrainedOptProblem(_constrainedOptProblem) {}

    double operator() (double x, int j) override {
        int numberConstraints = constrainedOptProblem->GetConstraintsNumber();
        if (j >= 1 && j <= numberConstraints) {
            return constrainedOptProblem->ComputeConstraint(j - 1, vector<double>{ x });
        } else if (j == numberConstraints + 1) {
            return constrainedOptProblem->ComputeFunction(vector<double>{ x });
        } else {
            return numeric_limits<double>::quiet_NaN();
        }
    }

    double operator() (vector<double> x, int j) override {
        int numberConstraints = constrainedOptProblem->GetConstraintsNumber();
        if (j >= 1 && j <= numberConstraints) {
            return constrainedOptProblem->ComputeConstraint(j - 1, x);
        } else if (j == numberConstraints + 1) {
            return constrainedOptProblem->ComputeFunction(x);
        } else {
            return numeric_limits<double>::quiet_NaN();
        }
    }
};

class FunctorFamily : public Functor {
public:
    IOptProblemFamily *optProblemFamily;
    int currentFunction;

    FunctorFamily(IOptProblemFamily *_optProblemFamily = nullptr) 
                  : Functor(), optProblemFamily(_optProblemFamily), currentFunction(0) {};

    double operator() (double x, int j) override {
        switch (j) {
            case 1: return (*optProblemFamily)[currentFunction]->ComputeFunction(vector<double>{ x });
            default: return numeric_limits<double>::quiet_NaN();
        }
    }

    double operator() (vector<double> x, int j) override {
        switch (j) {
            case 1: return (*optProblemFamily)[currentFunction]->ComputeFunction(x);
            default: return numeric_limits<double>::quiet_NaN();
        }
    }
};

class FunctorFamilyConstrained : public Functor {
public:
    IConstrainedOptProblemFamily *constrainedOptProblemFamily;
    int currentFunction;

    FunctorFamilyConstrained(IConstrainedOptProblemFamily *_constrainedOptProblemFamily = nullptr) : Functor(),
                             constrainedOptProblemFamily(_constrainedOptProblemFamily), currentFunction(0) {}

    double operator() (double x, int j) {
        int numberConstraints = (*constrainedOptProblemFamily)[currentFunction]->GetConstraintsNumber();
        if (j >= 1 && j <= numberConstraints) {
            return (*constrainedOptProblemFamily)[currentFunction]->ComputeConstraint(j - 1, vector<double>{ x });
        } else if (j == numberConstraints + 1) {
            return (*constrainedOptProblemFamily)[currentFunction]->ComputeFunction(vector<double>{ x });
        } else {
            return numeric_limits<double>::quiet_NaN();
        }
    }

    double operator() (vector<double> x, int j) {
        int numberConstraints = (*constrainedOptProblemFamily)[currentFunction]->GetConstraintsNumber();
        if (j >= 1 && j <= numberConstraints) {
            return (*constrainedOptProblemFamily)[currentFunction]->ComputeConstraint(j - 1, x);
        } else if (j == numberConstraints + 1) {
            return (*constrainedOptProblemFamily)[currentFunction]->ComputeFunction(x);
        } else {
            return numeric_limits<double>::quiet_NaN();
        }
    }
};

struct DataDirect {
    int numberFevals;
    vector<vector<double>> points;

    DataDirect() : numberFevals(0), points(0) {};
};

struct DataDirectOperationalCharacteristics : public DataDirect {
    Functor *functor;
    TypeConstraints type;
    bool converge;
    int minNumberFevals;
    double eps;

    DataDirectOperationalCharacteristics(Functor *_functor = nullptr, TypeConstraints _type = TypeConstraints::NoConstraints, double _eps = 0.01)
        : DataDirect(), functor(_functor), type(_type), converge(false), minNumberFevals(-1), eps(_eps) {}; 
};

struct TaskDirect {
    string name;

    direct_objective_func f;
    void *fData;
    int n;
    vector<double> A, B, XOpt;

    int maxFevals, maxIters;
    double magicEps;
    double volumeReltol, sigmaReltol;

    FILE *logfile;
    direct_algorithm algorithm;

    bool used;

    TaskDirect(string _name, direct_objective_func _f, void *_fData, int _n, vector<double> _A, vector<double> _B,
               vector<double> _XOpt, int _maxFevals, int _maxIters, double _magicEps, double _volumeReltol, double _sigmaReltol,
               FILE *_logfile, direct_algorithm _algorithm, bool _used = true) : name(_name), f(_f), fData(_fData), n(_n), A(_A),
               B(_B), XOpt(_XOpt), maxFevals(_maxFevals), maxIters(_maxIters), magicEps(_magicEps), volumeReltol(_volumeReltol),
               sigmaReltol(_sigmaReltol), logfile(_logfile), algorithm(_algorithm), used(_used) {};
}; 
*/
