#ifndef I_GENERAL_NUMERICAL_OPT_METHOD_H_
#define I_GENERAL_NUMERICAL_OPT_METHOD_H_

#include <vector>

#include <base_classes/opt_methods/IGeneralOptMethod.h>

using std::vector;

namespace opt {
    template <typename TrialType, typename PointType, typename OptProblemType,
              typename ParametersMethodType, typename ResultMethodType>
    class IGeneralNumericalOptMethod:
        public IGeneralOptMethod<OptProblemType, ParametersMethodType, ResultMethodType> {
    protected:
        vector<TrialType> trialPoints;

        int numberTrials, numberFevals;

        virtual TrialType newTrial(const PointType &x) = 0;
        virtual PointType selectNewPoint() = 0;

        virtual double estimateSolution(PointType &x) const = 0;

        virtual bool stopConditions() = 0;
        virtual bool stopConditionsTest() = 0;

    public:
        IGeneralNumericalOptMethod(const OptProblemType &_problem, const ParametersMethodType &_parameters):
            IGeneralOptMethod<OptProblemType, ParametersMethodType, ResultMethodType>(_problem, _parameters),
            trialPoints(0),
            numberTrials(0),
            numberFevals(0)
        {};

        void setAccuracy(double _accuracy) { parameters.accuracy = _accuracy; };
        double getAccuracy() const { return parameters.accuracy; };

        void setError(double _error) { parameters.error = _error; };
        double getError() const { return parameters.error; };

        void setMaxTrials(int _maxTrials) { parameters.maxTrials = _maxTrials; };
        int getMaxTrials() const { return parameters.maxTrials; };

        void setMaxFevals(int _maxFevals) { parameters.maxFevals = _maxFevals; };
        int getMaxFevals() const { return parameters.maxFevals; };

        void getTrialPoints(vector<TrialType> &_trialPoints) const { _trialPoints = trialPoints; };
        int getNumberTrialPoints() const { return trialPoints.size(); };

        using IGeneralOptMethod<OptProblemType, ParametersMethodType, ResultMethodType>::parameters;
    };
}

#endif // I_GENERAL_NUMERICAL_OPT_METHOD_H_
