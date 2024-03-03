#ifndef SOLVER_H_
#define SOLVER_H_

#include <iostream>
#include <algorithm>
#include <utility>
#include <vector>
#include <string>
#include <memory>

#include <general/classes/opt_methods/IGeneralOptMethod.h>
#include <general/classes/opt_methods/IGeneralNumericalOptMethod.h>
#include <gnuplot/OutputFile.h>
#include <MyMath.h>
#include <omp.h>

template <typename OptProblemType>
class Solver {
public:
    using Result = typename opt::IGeneralOptMethod<OptProblemType>::Result;
    using Report = typename opt::IGeneralOptMethod<OptProblemType>::IReport;

    Solver() = default;

    void solveTask(
        opt::IGeneralOptMethod<OptProblemType> &optMethod,
        const opt::Task<OptProblemType> &task, Result &result)
    {
        if (task.use) {
            optMethod.setProblem(task.problem);
            optMethod.setParameters(task.parameters);

            std::unique_ptr<Report> report(optMethod.createReport());

            double startTime = omp_get_wtime();
            optMethod.solve(result);
            double endTime = omp_get_wtime();
            double workTime = endTime - startTime;

            report->print(std::cout, task, result, workTime);
        }
    }

    void solveTasks(
        opt::IGeneralOptMethod<OptProblemType> &optMethod,
        const std::vector<opt::Task<OptProblemType>> &tasks)
    {
        size_t numberTasks = tasks.size();
        std::unique_ptr<Result> result(optMethod.createResult());

        double totalStartTime = omp_get_wtime();
        for (size_t i = 0; i < numberTasks; ++i) {
            result = solveTask(optMethod, tasks[i], *result);
            std::cout << "\n\n";
        }
        double totalEndTime = omp_get_wtime();

        std::cout << "Total time: " << totalEndTime - totalStartTime << "\n";
    }

    void solveTestTask(
        opt::IGeneralOptMethod<OptProblemType> &optMethod,
        const opt::Task<OptProblemType> &task, Result *result) {
        if (task.use) {
            optMethod.setProblem(task.problem);
            optMethod.setParameters(task.parameters);

            std::unique_ptr<Report> report(optMethod.createReport());

            double startTime = omp_get_wtime();
            optMethod.solveTest(result);
            double endTime = omp_get_wtime();
            double workTime = endTime - startTime;

            report->print(std::cout, task, result, workTime);
        }
    }

    void solveTestTasks(
        opt::IGeneralOptMethod<OptProblemType> &optMethod,
        const std::vector<opt::Task<OptProblemType>> &tasks)
    {
        size_t numberTasks = tasks.size();
        std::unique_ptr<Result> result(optMethod.createResult());

        double totalStartTime = omp_get_wtime();
        for (size_t i = 0; i < numberTasks; ++i) {
            solveTestTask(optMethod, tasks[i], result);
            std::cout << "\n\n";
        }
        double totalEndTime = omp_get_wtime();

        std::cout << "Total time: " << totalEndTime - totalStartTime << "\n";
    }

    template <typename TrialType>
    void calcOperationalCharacteristics(
        opt::IGeneralNumericalOptMethod<TrialType, OptProblemType> &optMethod,
        const opt::Task<OptProblemType> &task,
        size_t kStart, size_t kFinish, size_t kStep,
        std::vector<std::pair<size_t, double>> &operationalCharacteristics, double &workTime)
    {
        using Result = typename opt::IGeneralNumericalOptMethod<TrialType, OptProblemType>::Result;

        size_t numberPoints = (kFinish - kStart) / kStep + 1;
        operationalCharacteristics.resize(numberPoints);

        size_t familySize = task.problem.getFamilySize();
        std::unique_ptr<Result> result(static_cast<Result*>(optMethod.createResult()));
        std::vector<size_t> numberTrials(familySize);

        optMethod.setParameters(task.parameters);

        double startTime = omp_get_wtime();
        for (size_t i = 0; i < familySize; ++i) {
            task.problem.setProblemNumber(i);
            optMethod.setProblem(task.problem);

            numberTrials[i] = optMethod.solveTest(*result) ? result->numberTrials : kFinish + 1;
        }
        double endTime = omp_get_wtime();
        workTime = endTime - startTime;

        size_t numberSuccessful, k;
        for (size_t i = 0, k = kStart; k <= kFinish; ++i, k += kStep) {
            numberSuccessful = std::count_if(numberTrials.begin(), numberTrials.end(),
                                             [k] (double elem) { return elem <= k; });
            operationalCharacteristics[i] = std::pair<size_t, double>(k, (double)numberSuccessful / familySize);
        }
    }
};

#endif // SOLVER_H_
