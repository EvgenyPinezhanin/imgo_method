#ifndef FAMILY_PROBLEM_H_
#define FAMILY_PROBLEM_H_

#include <vector>
#include <string>
#include <limits>

#include <IConstrainedOptProblemFamily.hpp>

struct Problem {
    string name;
    string shortName;

    bool used;

    Problem(string _name = "", string _shortName = "", bool _used = true) 
            : name(_name), shortName(_shortName), used(_used) {};
};

struct FamilyProblem : public Problem {
    IConstrainedOptProblemFamily *optProblemFamily;

    FamilyProblem(string _name = "", IConstrainedOptProblemFamily *_optProblemFamily = nullptr,
                  string _shortName = "", bool _used = true) 
        : Problem(_name, _shortName, _used), optProblemFamily(_optProblemFamily) {};
};

/* struct ProblemSingle : public Problem {
    IGeneralOptProblem *optProblem;

    ProblemSingle(string _name, IGeneralOptProblem *_optProblem, TypeConstraints _type, string _shortName = "", bool _used = true) 
                  : Problem(_name, _type, _shortName, used), optProblem(_optProblem) {};
}; */

class Functor {
public:
    Functor() {};

    virtual double operator() (double x, int j) = 0;
    virtual double operator() (vector<double> x, int j) = 0; 
};

/* class FunctorSingle : public Functor {
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
}; */

/* class FunctorSingleConstrained : public Functor {
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
}; */

/* class FunctorFamily : public Functor {
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
}; */

class FunctorFamilyConstrained : public Functor {
public:
    IConstrainedOptProblemFamily *constrainedOptProblemFamily;
    int currentFunction;

    FunctorFamilyConstrained(IConstrainedOptProblemFamily *_constrainedOptProblemFamily = nullptr) : Functor(),
                             constrainedOptProblemFamily(_constrainedOptProblemFamily), currentFunction(0) {}

    double operator() (double x, int j) {
        int numberConstraints = (*constrainedOptProblemFamily)[currentFunction]->GetConstraintsNumber();
        if (j >= 0 && j < numberConstraints) {
            return (*constrainedOptProblemFamily)[currentFunction]->ComputeConstraint(j, vector<double>{ x });
        } else if (j == numberConstraints) {
            return (*constrainedOptProblemFamily)[currentFunction]->ComputeFunction(vector<double>{ x });
        } else {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }

    double operator() (vector<double> x, int j) override {
        int numberConstraints = (*constrainedOptProblemFamily)[currentFunction]->GetConstraintsNumber();
        if (j >= 0 && j < numberConstraints) {
            return (*constrainedOptProblemFamily)[currentFunction]->ComputeConstraint(j, x);
        } else if (j == numberConstraints) {
            return (*constrainedOptProblemFamily)[currentFunction]->ComputeFunction(x);
        } else {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }
};

/*
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


#endif // FAMILY_PROBLEM_H_
