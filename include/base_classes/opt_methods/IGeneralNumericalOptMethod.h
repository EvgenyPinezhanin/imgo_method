#ifndef I_GENERAL_NUMERICAL_OPT_METHOD_H_
#define I_GENERAL_NUMERICAL_OPT_METHOD_H_

#include <vector>

#include <base_classes/opt_methods/IGeneralOptMethod.h>
#include <base_structures/GeneralParametersNumericalOptMethod.h>

namespace opt {
    template <typename TrialType, typename PointType, typename OptProblemType, typename ResultMethodType>
    class IGeneralNumericalOptMethod:
        public IGeneralOptMethod<OptProblemType, GeneralParametersNumericalOptMethod, ResultMethodType> {
    protected:
        std::vector<TrialType> trialPoints;

        double accuracy, error;
        int numberTrials, maxTrials;
        int numberFevals, maxFevals;

        virtual TrialType newTrial(const PointType &x) = 0;
        virtual PointType selectNewPoint() = 0;

        virtual double estimateSolution(PointType &x) const = 0;

        virtual bool stopConditions() = 0;
        virtual bool stopConditionsTest() = 0;

    public:
        IGeneralNumericalOptMethod(const OptProblemType &_problem,
                                   const GeneralParametersNumericalOptMethod &parameters):
            IGeneralOptMethod<OptProblemType, GeneralParametersNumericalOptMethod, ResultMethodType>(_problem),
            trialPoints(),
            accuracy(parameters.accuracy),
            error(parameters.error),
            numberTrials(0),
            maxTrials(parameters.maxTrials),
            numberFevals(0),
            maxFevals(parameters.maxFevals)
        {};

        void setParameters(const GeneralParametersNumericalOptMethod &parameters) override {
            accuracy = parameters.accuracy;
            error = parameters.error;
            maxTrials = parameters.maxTrials;
            maxFevals = parameters.maxFevals;
        }
        void getParameters(GeneralParametersNumericalOptMethod &parameters) const override {
            parameters = GeneralParametersNumericalOptMethod(accuracy, error, numberTrials, numberFevals);
        }

        void setAccuracy(double _accuracy) { accuracy = _accuracy; };
        double getAccuracy() const { return accuracy; };

        void setError(double _error) { error = _error; };
        double getError() const { return error; };

        void setMaxTrials(int _maxTrials) { maxTrials = _maxTrials; };
        int getMaxTrials() const { return maxTrials; };

        void setMaxFevals(int _maxFevals) { maxFevals = _maxFevals; };
        int getMaxFevals() const { return maxFevals; };

        void getTrialPoints(std::vector<TrialType> &_trialPoints) const { _trialPoints = trialPoints; };
        int getNumberTrialPoints() const { return trialPoints.size(); };
    };
}

#endif // I_GENERAL_NUMERICAL_OPT_METHOD_H_
