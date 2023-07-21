#ifndef I_GENERAL_CONSTRAINED_OPT_TASK_H_
#define I_GENERAL_CONSTRAINED_OPT_TASK_H_

#include <base_classes/opt_problems/IGeneralOptProblem.h>

namespace om {
    template <typename ObjectiveFunctionType, typename SearchAreaType, typename PointType>
    class IGeneralConstrainedOptProblem : public IGeneralOptProblem<ObjectiveFunctionType, SearchAreaType, PointType> {
    protected:
        int numberConstraints;

    public:
        IGeneralConstrainedOptProblem(const ObjectiveFunctionType &_objFunction, int _numberConstraints,
                                      const SearchAreaType &_area, const PointType &_optPoint)
                                     : IGeneralOptProblem<ObjectiveFunctionType, SearchAreaType, PointType>(_objFunction,
                                     _area, _optPoint), numberConstraints(_numberConstraints) {};

        void setNumberConstraints(double _numberConstraints) { numberConstraints = _numberConstraints; };
        int getNumberConstraints() const { return numberConstraints; };

        virtual double computeConstraints(const PointType &x, int index) const = 0;
    };
}

#endif // I_GENERAL_CONSTRAINED_OPT_TASK_H_
