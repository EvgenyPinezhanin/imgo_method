#ifndef SOLVER_H_
#define SOLVER_H_

#include <iostream>
#include <vector>
#include <string>

#include <base_classes/opt_methods/IGeneralOptMethod.h>
#include <gnuplot/OutputFile.h>
#include <omp.h>

template <typename OptMethodType>
class Solver {
protected:
    OptMethodType optMethod;
    typename OptMethodType::Parameters parameters;
    typename OptMethodType::Result result;
    typename OptMethodType::Report report;

    OutputFile trialsFile;

    std::vector<typename OptMethodType::Trial> trials;
    std::vector<typename OptMethodType::OptProblem::Point> optimalPoints;
    typename OptMethodType::OptProblem::Point point;

public:
    Solver(const OptMethodType &_optMethod):
        optMethod(_optMethod),
        result(),
        report(),
        trialsFile()
    {};

    void solveTask(const opt::Task<typename OptMethodType::OptProblem> &task, const std::string &saveDirectory) {
        if (task.use) {
            std::string blockName;
            task.problem.getBlockName(blockName);
            trialsFile.open(saveDirectory + blockName + "_" + std::to_string(task.problem.getProblemNumber()));
            if (!trialsFile.isOpen()) std::cerr << "trials_file opening error\n";

            optMethod.setProblem(task.problem);
            optMethod.setParameters(task.parameters);

            double startTime = omp_get_wtime();
            optMethod.solve(result);
            double endTime = omp_get_wtime();
            double workTime = endTime - startTime;

            report.print(std::cout, task, result, workTime);

            task.problem.getOptimalPoints(optimalPoints);
            trialsFile.addPoints(optimalPoints, task.problem.getOptimalValue());

            trialsFile.addPoint(result.point, result.value);

            optMethod.getTrialPoints(trials);
            trialsFile.addPoints(trials);

            trialsFile.close();
        }
    }

    void solveTasks(const std::vector<opt::Task<typename OptMethodType::OptProblem>> &tasks, const std::string &saveDirectory) {
        size_t numberTasks = tasks.size();

        double totalStartTime = omp_get_wtime();
        for (size_t i = 0; i < numberTasks; ++i) {
            solveTask(tasks[i], saveDirectory);
            std::cout << "\n\n";
        }
        double totalEndTime = omp_get_wtime();

        std::cout << "Total time: " << totalEndTime - totalStartTime << "\n";
    }

    void solveTestTask(const opt::Task<typename OptMethodType::OptProblem> &task, const std::string &saveDirectory) {
        if (task.use) {
            trialsFile.open(saveDirectory + task.blockName + "_" + std::to_string(task.functionNumber));
            if (!trialsFile.isOpen()) std::cerr << "trials_file opening error\n";

            optMethod.setProblem(task.problem);
            optMethod.setParameters(task.parameters);

            double startTime = omp_get_wtime();
            optMethod.solveTest(result);
            double endTime = omp_get_wtime();
            double workTime = endTime - startTime;

            report.print(std::cout, task, result, workTime);

            task.problem.getOptimalPoints(optimalPoints);
            trialsFile.addPoints(optimalPoints, task.problem.getOptimalValue());

            trialsFile.addPoint(result.point, result.value);

            optMethod.getTrialPoints(trials);
            trialsFile.addPoints(trials);

            trialsFile.close();
        }
    }

    void solveTestTasks(const std::vector<opt::Task<typename OptMethodType::OptProblem>> &tasks, const std::string &saveDirectory) {
        size_t numberTasks = tasks.size();

        double totalStartTime = omp_get_wtime();
        for (size_t i = 0; i < numberTasks; ++i) {
            solveTestTask(tasks[i], saveDirectory);
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
