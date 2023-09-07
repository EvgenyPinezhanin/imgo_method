#ifndef I_GENERAL_REPORT_H_
#define I_GENERAL_REPORT_H_

#include <ostream>
#include <iomanip>

#include <base_structures/GeneralParametersNumericalOptMethod.h>

namespace opt {
    template <typename TaskType, typename OptProblemType, typename ResultMethodType>
    class IGeneralReport {
    protected:
        virtual void printOptProblem(std::ostream &stream, const OptProblemType &optProblem) const = 0;
        virtual void printMethodParameters(std::ostream &stream,
                                           const GeneralParametersNumericalOptMethod &parameters) const = 0;
        virtual void printResultMethod(std::ostream &stream, const ResultMethodType &result) const = 0;
        virtual void printErrorEstimate(std::ostream &stream, const OptProblemType &optProblem,
                                        const ResultMethodType &result) const = 0;

    public:
        IGeneralReport() {};

        void print(std::ostream &stream, const TaskType &task, const ResultMethodType &result, double workTime) const {
            const auto defaultPrecision = stream.precision();
            stream << std::setprecision(10);

            stream << "Task: " << task.name << "\n";
            printOptProblem(stream, task.optProblem);
            printMethodParameters(stream, task.parameters);
            printResultMethod(stream, result);
            printErrorEstimate(stream, task.optProblem, result);
            stream << "Time: " << workTime << std::endl;

            stream << std::setprecision(defaultPrecision);
        }
    };
}

#endif // I_GENERAL_REPORT_H_
