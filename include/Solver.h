#ifndef SOLVER_H_
#define SOLVER_H_

#include <iostream>
#include <vector>
#include <string>
#include <memory>

#include <base_classes/opt_methods/IGeneralOptMethod.h>
#include <gnuplot/OutputFile.h>
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

    void calcOperationalCharacteristics() {
/*         cout << "r = " << r[i][j] << endl;
            imgo.setR(r[i][j]);

            double startTime = omp_get_wtime();
            for (int k = 0; k < numberFunctions; k++) {
                functor.currentFunction = k;
                imgo.setF(function<double(double, int)>(functor));
                imgo.setMaxTrials(Kmax);
                (*functor.optProblemFamily)[k]->GetBounds(A, B);
                imgo.setAB(A[0], B[0]);

                if (imgo.solveTest((*functor.optProblemFamily)[k]->GetOptimumPoint()[0], numberTrials, numberFevals)) {
                    numberTrialsArray[i][k] = numberTrials;
                } else {
                    numberTrialsArray[i][k] = Kmax + 1;
                }
            }
            for (int k = K0; k <= Kmax; k += Kstep) {
                numberSuccessful = (int)count_if(numberTrialsArray[i].begin(), numberTrialsArray[i].end(),
                                                 [k] (double elem) { return elem <= k; });
                cout << "K = " << k << " success rate = " << (double)numberSuccessful / numberFunctions << endl;
                ofstr << k << " " << (double)numberSuccessful / numberFunctions << endl;
            }
            ofstr << endl << endl;
            double endTime = omp_get_wtime();
            double workTime = endTime - startTime;
            cout << "Time: " << workTime << endl; */
    }
};

#endif // SOLVER_H_
