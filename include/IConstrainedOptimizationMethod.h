#ifndef I_CONSTRAINED_OBJECTIVE_FUNCTION_H_
#define I_CONSTRAINED_OBJECTIVE_FUNCTION_H_

#include <vector>

using std::vector;

class IConstrainedOptimizationMethod {
protected:
    int numberConstraints;

    virtual double computeConstraint(vector<double> X, int index) const = 0;

public:
    IConstrainedOptimizationMethod(int _numberConstraints) : numberConstraints(_numberConstraints) {};

    void setNumberConstraints(int _numberConstraints) { numberConstraints = _numberConstraints; };
    int getNumberConstraints() const { return numberConstraints; };
};

#endif // I_CONSTRAINED_OBJECTIVE_FUNCTION_H_
