#ifndef _I_GENERAL_NUMERICAL_OPT_METHOD_H_
#define _I_GENERAL_NUMERICAL_OPT_METHOD_H_

#include <vector>

#include <general/classes/opt_methods/IGeneralOptMethod.h>

namespace opt {
    static size_t initTrialPointsSize = 100;

    template <typename TrialType, typename OptProblemType>
    class IGeneralNumericalOptMethod : public IGeneralOptMethod<OptProblemType> {
    public:
        using GeneralMethod = IGeneralOptMethod<OptProblemType>;
        using Trial = TrialType;

        struct Parameters : public GeneralMethod::Parameters {
            double accuracy, error;
            size_t maxTrials, maxFevals;

            Parameters(double _accuracy = 0.001, double _error = 0.001, size_t _maxTrials = 1000, size_t _maxFevals = 1000):
                accuracy(_accuracy), error(_error), maxTrials(_maxTrials), maxFevals(_maxFevals) {};
        };

        using StoppingCondition = size_t;

        struct StoppingConditions {
            static const StoppingCondition ACCURACY = 0;
            static const StoppingCondition ERROR = 1;
            static const StoppingCondition MAXTRIALS = 2;
            static const StoppingCondition MAXFEVALS = 3;
        };

        struct Result : public GeneralMethod::Result {
            size_t numberTrials, numberFevals;
            double resultingAccuracy;

            StoppingCondition stoppingCondition;

            Result(const typename OptProblemType::Point &_point = typename OptProblemType::Point(),
                   double _value = 0.0, size_t _numberTrials = 0, size_t _numberFevals = 0, double _resultingAccuracy = 0.0,
                   StoppingCondition _stoppingCondition = StoppingConditions::ACCURACY):
                GeneralMethod::Result(_point, _value), numberTrials(_numberTrials), numberFevals(_numberFevals),
                resultingAccuracy(_resultingAccuracy), stoppingCondition(_stoppingCondition) {};
            ~Result() override {};
        };

        class IReport : public GeneralMethod::IReport {
        protected:
            // TODO: add function for printing bounds
            void printMethodParameters(std::ostream &stream,
                                       const typename GeneralMethod::Parameters &parameters) const override;
            void printResultMethod(std::ostream &stream,
                                   const typename GeneralMethod::Result &result) const override;

        public:
            IReport() : GeneralMethod::IReport() {};
            ~IReport() override {};

            virtual void printStopCondition(std::ostream &stream, StoppingCondition stoppingCondition) const;
        };

    protected:
        std::vector<TrialType> trialPoints;

        double accuracy, resultingAccuracy, error;
        size_t numberTrials, maxTrials;
        size_t numberFevals, maxFevals;

        StoppingCondition stoppingCondition;

        static constexpr double epsilon = 1e-14;

        virtual TrialType newTrial(const typename OptProblemType::Point &x) = 0;
        virtual typename OptProblemType::Point selectNewPoint() = 0;

        virtual double estimateSolution(typename OptProblemType::Point &x) const = 0;

        virtual bool stopConditions() = 0;
        virtual bool stopConditionsTest() = 0;

        void setResult(typename GeneralMethod::Result &result) const override;

    public:
        IGeneralNumericalOptMethod(const OptProblemType &_problem, const Parameters &parameters)
            : GeneralMethod(_problem), trialPoints(initTrialPointsSize), accuracy(parameters.accuracy),
              resultingAccuracy(0.0), error(parameters.error), numberTrials(0), maxTrials(parameters.maxTrials),
              numberFevals(0), maxFevals(parameters.maxFevals), stoppingCondition(StoppingConditions::ACCURACY) {};

        void setParameters(const typename GeneralMethod::Parameters &parameters) override {
            auto parametersCast = static_cast<const Parameters&>(parameters);
            accuracy = parametersCast.accuracy;
            error = parametersCast.error;
            maxTrials = parametersCast.maxTrials;
            maxFevals = parametersCast.maxFevals;
        }
        void getParameters(typename GeneralMethod::Parameters &parameters) const override {
            parameters = Parameters(accuracy, error, maxTrials, maxFevals);
        }

        void setAccuracy(double _accuracy) { accuracy = _accuracy; };
        double getAccuracy() const { return accuracy; };

        void setError(double _error) { error = _error; };
        double getError() const { return error; };

        void setMaxTrials(size_t _maxTrials) { maxTrials = _maxTrials; };
        size_t getMaxTrials() const { return maxTrials; };

        void setMaxFevals(size_t _maxFevals) { maxFevals = _maxFevals; };
        size_t getMaxFevals() const { return maxFevals; };

        void getTrialPoints(std::vector<TrialType> &_trialPoints) const { _trialPoints = trialPoints; };
        size_t getNumberTrialPoints() const { return trialPoints.size(); };

        typename GeneralMethod::Result* createResult() const override { return new Result(); };
        using GeneralMethod::createReport;
    };

    template <typename TrialType, typename OptProblemType>
    void IGeneralNumericalOptMethod<TrialType, OptProblemType>::IReport::printMethodParameters(
        std::ostream &stream, const typename GeneralMethod::Parameters &parameters) const
    {
        auto parametersCast = static_cast<const Parameters&>(parameters);
        stream << "Method Parameters:" << "\n";
        stream << "Maximum of trials = " << parametersCast.maxTrials << "\n";
        stream << "Maximum of fevals = " << parametersCast.maxFevals << "\n";
        stream << "Accuracy = " << parametersCast.accuracy << "\n";
    }

    template <typename TrialType, typename OptProblemType>
    void IGeneralNumericalOptMethod<TrialType, OptProblemType>::IReport::printResultMethod(
        std::ostream &stream, const typename GeneralMethod::Result &result) const
    {
        auto resultCast = static_cast<const Result&>(result);
        stream << "Result of method:" << "\n";
        stream << "Number of trials = " << resultCast.numberTrials << "\n";
        stream << "Number of fevals = " << resultCast.numberFevals << "\n";
        stream << "Resulting accuracy = " << resultCast.resultingAccuracy << "\n";

        stream << "Stop condition: ";
        this->printStopCondition(stream, resultCast.stoppingCondition);
        stream << "\n";

        stream << "X = ";
        this->printPoint(stream, resultCast.point);
        stream << "\n";
        stream << "f(X) = " << resultCast.value << "\n";
    }

    template <typename TrialType, typename OptProblemType>
    void IGeneralNumericalOptMethod<TrialType, OptProblemType>::IReport::printStopCondition(
        std::ostream &stream, StoppingCondition stoppingCondition) const
    {
        switch (stoppingCondition) {
            case StoppingConditions::ACCURACY:
                stream << "Accuracy";
                break;
            case StoppingConditions::ERROR:
                stream << "Error";
                break;
            case StoppingConditions::MAXTRIALS:
                stream << "Max trials";
                break;
            case StoppingConditions::MAXFEVALS:
                stream << "Max fevals";
                break;
            default: break;
        }
    }


    template <typename TrialType, typename OptProblemType>
    void IGeneralNumericalOptMethod<TrialType, OptProblemType>::setResult(
        typename GeneralMethod::Result &result) const
    {
        typename OptProblemType::Point x;

        result.value = estimateSolution(x);
        result.point = x;

        auto& resultCast = static_cast<Result&>(result);
        resultCast.numberTrials = numberTrials;
        resultCast.numberFevals = numberFevals;
        resultCast.resultingAccuracy = resultingAccuracy;
        resultCast.stoppingCondition = stoppingCondition;
    }
}

#endif // _I_GENERAL_NUMERICAL_OPT_METHOD_H_
