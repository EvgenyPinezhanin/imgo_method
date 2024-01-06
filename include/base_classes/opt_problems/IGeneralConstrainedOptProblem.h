#ifndef I_GENERAL_CONSTRAINED_OPT_TASK_H_
#define I_GENERAL_CONSTRAINED_OPT_TASK_H_

#include <base_classes/opt_problems/IGeneralOptProblem.h>

namespace opt {
    template <typename ObjectiveFunctionType, typename SearchAreaType, typename PointType>
    class IGeneralConstrainedOptProblem : public IGeneralOptProblem<ObjectiveFunctionType, SearchAreaType, PointType> {
    protected:
        size_t numberConstraints;

    public:
        IGeneralConstrainedOptProblem(const ObjectiveFunctionType &_objFunction, const SearchAreaType &_area,
                                      size_t _numberConstraints, const std::vector<PointType> &_optPoint,  double _optimalValue)
            : IGeneralOptProblem<ObjectiveFunctionType, SearchAreaType, PointType>(_objFunction,
              _area, _optPoint, _optimalValue), numberConstraints(_numberConstraints) {};

        virtual size_t getNumberConstraints() const { return numberConstraints; };

        virtual double computeConstraint(const PointType &x, int index) const = 0;
    };
}

#endif // I_GENERAL_CONSTRAINED_OPT_TASK_H_
