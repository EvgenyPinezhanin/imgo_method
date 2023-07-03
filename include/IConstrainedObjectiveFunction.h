#ifndef I_CONSTRAINED_OBJECTIVE_FUNCTION_H_
#define I_CONSTRAINED_OBJECTIVE_FUNCTION_H_

#include <vector>

#include <IObjectiveFunction.h>

using std::vector;

template <typename ObjectiveFunction>
class IConstrainedObjectiveFunction : public IObjectiveFunction<ObjectiveFunction> {
protected:
    int numberConstraints;

    virtual double computeConstraint(vector<double> X, int index) const = 0;

public:
    IConstrainedObjectiveFunction(ObjectiveFunction _objFunction, int _numberConstraints)
                                  : objFunction(_objFunction), numberConstraints(_numberConstraints) {};

    void setNumberConstraints(int _numberConstraints) { numberConstraints = _numberConstraints; };

    int getNumberConstraints() const { return numberConstraints; };
};

#endif // I_CONSTRAINED_OBJECTIVE_FUNCTION_H_
