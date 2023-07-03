#ifndef I_OBJECTIVE_FUNCTION_H_
#define I_OBJECTIVE_FUNCTION_H_

#include <vector>

using std::vector;

template <typename ObjectiveFunction>
class IObjectiveFunction {
protected:
    ObjectiveFunction objFunction;

    virtual double compute(vector<double> X) const = 0;

public:
    IObjectiveFunction(ObjectiveFunction _objFunction) : objFunction(_objFunction) {};

    void setF(ObjectiveFunction _objFunction) { objFunction = _objFunction; };

    ObjectiveFunction getF() const { return objFunction; };
};

#endif // I_OBJECTIVE_FUNCTION_H_
