#ifndef _I_CONSTRAINED_OPT_PROBLEM_H_
#define _I_CONSTRAINED_OPT_PROBLEM_H_

#include <vector>

namespace opt {
    template <typename PointType>
    class IConstrainedOptProblem {
    protected:
        size_t numberConstraints;

    public:
        IConstrainedOptProblem(size_t _numberConstraints) : numberConstraints(_numberConstraints) {};

        virtual size_t getNumberConstraints() const { return numberConstraints; };
        virtual double computeConstraintFunction(const PointType &x, size_t index) const = 0;
    };
}

#endif // _I_CONSTRAINED_OPT_PROBLEM_H_
