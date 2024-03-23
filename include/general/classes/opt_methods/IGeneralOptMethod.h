#ifndef _I_GENERAL_OPT_METHOD_H_
#define _I_GENERAL_OPT_METHOD_H_

#include <iostream>
#include <iomanip>

namespace opt {
    template <typename OptProblemType>
    class IGeneralOptMethod;

    // TODO: make Parameters as template parameter
    template <typename OptProblemType>
    struct Task {
        using Parameters = typename IGeneralOptMethod<OptProblemType>::Parameters;

        std::string name;
        OptProblemType problem;
        Parameters &parameters;
        bool use;

        Task(const std::string &_name, const OptProblemType _problem,
             Parameters &_parameters, bool _use = true)
            : name(_name), problem(_problem), parameters(_parameters), use(_use) {};
    };

    template <typename OptProblemType>
    class IGeneralOptMethod {
    public:
        using OptProblem = OptProblemType;

        struct Parameters {
            Parameters() = default;
        };

        struct Result {
            typename OptProblemType::Point point;
            double value;

            Result(const typename OptProblemType::Point _point = typename OptProblemType::Point(),
                   double _value = 0.0)
                : point(_point), value(_value) {};
            virtual ~Result() {};
        };

        class IReport {
        protected:
            virtual void printOptProblem(std::ostream &stream, const OptProblemType &optProblem) const = 0;
            virtual void printMethodParameters(std::ostream &stream, const Parameters &parameters) const = 0;
            virtual void printResultMethod(std::ostream &stream, const Result &result) const = 0;
            virtual void printErrorEstimate(std::ostream &stream, const OptProblemType &optProblem,
                                            const Result &result) const = 0;

        public:
            IReport() = default;
            virtual ~IReport() {};

            virtual void printPoint(std::ostream &stream,
                                    const typename OptProblemType::Point &point) const = 0;
            virtual void printPoints(std::ostream &stream,
                                     const std::vector<typename OptProblemType::Point> &points) const;

            void print(std::ostream &stream, const Task<OptProblemType> &task, const Result &result, double workTime) const;
        };

    protected:
        OptProblemType problem;

        virtual void setResult(Result &result) const = 0;

    public:
        IGeneralOptMethod(const OptProblemType &_problem) : problem(_problem) {};

        void setProblem(const OptProblemType &_problem) { problem = _problem; };
        OptProblemType getProblem() const { return problem; };

        virtual void setParameters(const Parameters &parameters) = 0;
        virtual void getParameters(Parameters &parameters) const = 0;

        virtual Result* createResult() const { return new Result(); };
        virtual IReport* createReport() const = 0;

        virtual void solve(Result &result) = 0;
        virtual bool solveTest(Result &result) = 0;
    };

    template <typename OptProblemType>
    void IGeneralOptMethod<OptProblemType>::IReport::printPoints(
        std::ostream &stream, const std::vector<typename OptProblemType::Point> &points) const
    {
        stream << "[";
        printPoint(stream, points[0]);
        size_t numberPoints = points.size();
        for (size_t i = 1; i < numberPoints; ++i) {
            stream << "; ";
            printPoint(stream, points[i]);
        }
        stream << "]";
    }

    template <typename OptProblemType>
    void IGeneralOptMethod<OptProblemType>::IReport::print(
        std::ostream &stream, const Task<OptProblemType> &task, const Result &result, double workTime) const
    {
        const auto defaultPrecision = stream.precision();
        stream << std::setprecision(10);

        stream << "Task: " << task.name << "\n";
        printOptProblem(stream, task.problem);
        printMethodParameters(stream, task.parameters);
        printResultMethod(stream, result);
        printErrorEstimate(stream, task.problem, result);
        stream << "Time: " << workTime << std::endl;

        stream << std::setprecision(defaultPrecision);
    }
}

#endif // _I_GENERAL_OPT_METHOD_H_
