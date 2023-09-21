#ifndef I_GENERAL_NUMERICAL_OPT_METHOD_H_
#define I_GENERAL_NUMERICAL_OPT_METHOD_H_

#include <vector>

#include <base_classes/opt_methods/IGeneralOptMethod.h>

namespace opt {
    static size_t initTrialPointsSize = 100;

    template <typename TrialType, typename OptProblemType>
    class IGeneralNumericalOptMethod : public IGeneralOptMethod<OptProblemType> {
    public: 
        using GeneralMethod = IGeneralOptMethod<OptProblemType>;
        using Trial = TrialType;

        struct Parameters : public GeneralMethod::Parameters {
            double accuracy, error;
            int maxTrials, maxFevals;

            Parameters(double _accuracy = 0.001, double _error = 0.001, int _maxTrials = 1000, int _maxFevals = 1000):
                accuracy(_accuracy),
                error(_error),
                maxTrials(_maxTrials),
                maxFevals(_maxFevals)
            {};
        };

        using StoppingCondition = int;

        struct StoppingConditions {
            static const StoppingCondition ACCURACY = 0;
            static const StoppingCondition ERROR = 1;
            static const StoppingCondition MAXTRIALS = 2;
            static const StoppingCondition MAXFEVALS = 3;
        };

        struct Result : public GeneralMethod::Result {
            int numberTrials, numberFevals;

            StoppingCondition stoppingCondition;

            Result(const typename OptProblemType::Point &_point = typename OptProblemType::Point(),
                   double _value = 0.0, int _numberTrials = 0, int _numberFevals = 0,
                   StoppingCondition _stoppingCondition = StoppingConditions::ACCURACY):
                GeneralMethod::Result(_point, _value),
                numberTrials(_numberTrials),
                numberFevals(_numberFevals),
                stoppingCondition(_stoppingCondition)
            {};
        };

        class IReport : public GeneralMethod::IReport {
        protected:
            void printMethodParameters(std::ostream &stream,
                                       const typename GeneralMethod::Parameters &parameters) const override
            {
                auto parametersCast = static_cast<const Parameters&>(parameters);
                stream << "Method Parameters:" << "\n";
                stream << "Maximum of trials = " << parametersCast.maxTrials << "\n";
                stream << "Maximum of fevals = " << parametersCast.maxFevals << "\n";
                stream << "Accuracy = " << parametersCast.accuracy << "\n";
            }
            void printResultMethod(std::ostream &stream, const typename GeneralMethod::Result &result) const override {
                auto resultCast = static_cast<const Result&>(result);
                stream << "Result of method:" << "\n";
                stream << "Number of trials = " << resultCast.numberTrials << "\n";
                stream << "Number of fevals = " << resultCast.numberFevals << "\n";
                stream << "X = " << resultCast.point << "\n";
                stream << "f(X) = " << resultCast.value << "\n";
            }

        public:
            IReport():
                GeneralMethod::IReport()
            {};
        };

    protected:
        std::vector<TrialType> trialPoints;

        double accuracy, error;
        int numberTrials, maxTrials;
        int numberFevals, maxFevals;

        StoppingCondition stoppingCondition;
        
        virtual TrialType newTrial(const typename OptProblemType::Point &x) = 0;
        virtual typename OptProblemType::Point selectNewPoint() = 0;

        virtual double estimateSolution(typename OptProblemType::Point &x) const = 0;

        virtual bool stopConditions() = 0;
        virtual bool stopConditionsTest() = 0;

        void setResult(typename GeneralMethod::Result &result) const override {
            typename OptProblemType::Point x;

            result.value = estimateSolution(x);
            result.point = x;

            Result& resultNum = static_cast<Result&>(result);
            resultNum.numberTrials = numberTrials;
            resultNum.numberFevals = numberFevals;
            resultNum.stoppingCondition = stoppingCondition;
        }

    public:
        IGeneralNumericalOptMethod(const OptProblemType &_problem,
                                   const Parameters &parameters):
            GeneralMethod(_problem),
            trialPoints(initTrialPointsSize),
            accuracy(parameters.accuracy),
            error(parameters.error),
            numberTrials(0),
            maxTrials(parameters.maxTrials),
            numberFevals(0),
            maxFevals(parameters.maxFevals),
            stoppingCondition(StoppingConditions::ACCURACY)
        {};

        void setParameters(const typename GeneralMethod::Parameters &parameters) override {
            auto parametersTmp = static_cast<const Parameters&>(parameters);
            accuracy = parametersTmp.accuracy;
            error = parametersTmp.error;
            maxTrials = parametersTmp.maxTrials;
            maxFevals = parametersTmp.maxFevals;
        }
        void getParameters(typename GeneralMethod::Parameters &parameters) const override {
            parameters = Parameters(accuracy, error, maxTrials, maxFevals);
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
