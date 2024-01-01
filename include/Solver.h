#ifndef SOLVER_H_
#define SOLVER_H_

#include <iostream>
#include <vector>
#include <string>

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

    void solveTask(const typename OptMethodType::Task &task, const std::string &saveDirectory) {
        if (task.use) {
            trialsFile.open(saveDirectory + task.blockName + "_" + std::to_string(task.functionNumber));
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

    void solveTasks(const std::vector<typename OptMethodType::Task> &tasks, const std::string &saveDirectory) {
        size_t numberTasks = tasks.size();

        double totalStartTime = omp_get_wtime();
        for (size_t i = 0; i < numberTasks; ++i) {
            solveTask(tasks[i], saveDirectory);
            std::cout << "\n\n";
        }
        double totalEndTime = omp_get_wtime();

        std::cout << "Total time: " << totalEndTime - totalStartTime << "\n";
    }

    void solveTestTask(const typename OptMethodType::Task &task, const std::string &saveDirectory) {
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

    void solveTestTasks(const std::vector<typename OptMethodType::Task> &tasks, const std::string &saveDirectory) {
        size_t numberTasks = tasks.size();

        double totalStartTime = omp_get_wtime();
        for (size_t i = 0; i < numberTasks; ++i) {
            solveTestTask(tasks[i], saveDirectory);
            std::cout << "\n\n";
        }
        double totalEndTime = omp_get_wtime();

        std::cout << "Total time: " << totalEndTime - totalStartTime << "\n";
    }
};

#endif // SOLVER_H_
