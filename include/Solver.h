#ifndef SOLVER_H_
#define SOLVER_H_

#include <iostream>
#include <vector>
#include <string>

#include <test_opt_problems/onedimensional_opt_problems.h>
#include <gnuplot/output_file.h>
#include <omp.h>

template <typename OptMethodType>
class Solver {
protected:
    OptMethodType optMethod;
    typename OptMethodType::Parameters parameters;
    typename OptMethodType::Result result;
    typename OptMethodType::Report report;

    output_file trials_file;

    vector<typename OptMethodType::Trial> trials;
    vector<typename OptMethodType::OptProblem::Point> optimalPoints;
    typename OptMethodType::OptProblem::Point point;

public:
    Solver(const OptMethodType &_optMethod):
        optMethod(_optMethod),
        result(),
        report(),
        trials_file()
    {};

    void solveTask(const typename OptMethodType::Task &task, const string &saveDirectory) {
        if (task.use) {
            trials_file.open(saveDirectory + blockNames[task.blockNumber] + 
                             "_" + std::to_string(task.functionNumber));
            if (!trials_file.is_open()) std::cerr << "trials_file opening error\n";

            optMethod.setProblem(task.problem);
            optMethod.setParameters(task.parameters);

            double startTime = omp_get_wtime();
            optMethod.solve(result);
            double endTime = omp_get_wtime();
            double workTime = endTime - startTime;

            report.print(std::cout, task, result, workTime);

            task.problem.getOptimalPoints(optimalPoints);
            trials_file.add_points(optimalPoints, task.problem.getOptimalValue());

            trials_file.add_point(result.point, result.value);

            optMethod.getTrialPoints(trials);
            trials_file.add_points(trials);

            trials_file.close();
        }
    }

    void solveTasks(const vector<typename OptMethodType::Task> &tasks, const string &saveDirectory) {
        int numberTasks = (int)tasks.size();

        double totalStartTime = omp_get_wtime();
        for (int i = 0; i < numberTasks; i++) {
            solveTask(tasks[i], saveDirectory);
            std::cout << "\n\n";
        }
        double totalEndTime = omp_get_wtime();

        std::cout << "Total time: " << totalEndTime - totalStartTime << "\n";
    }

    void solveTestTask(const typename OptMethodType::Task &task, const string &saveDirectory) {
        if (task.use) {
            trials_file.open(saveDirectory + blockNames[task.blockNumber] + 
                             "_" + std::to_string(task.functionNumber));
            if (!trials_file.is_open()) std::cerr << "trials_file opening error\n";

            optMethod.setProblem(task.problem);
            optMethod.setParameters(task.parameters);

            double startTime = omp_get_wtime();
            optMethod.solveTest(result);
            double endTime = omp_get_wtime();
            double workTime = endTime - startTime;

            report.print(std::cout, task, result, workTime);

            task.problem.getOptimalPoints(optimalPoints);
            trials_file.add_points(optimalPoints, task.problem.getOptimalValue());

            trials_file.add_point(result.point, result.value);

            optMethod.getTrialPoints(trials);
            trials_file.add_points(trials);

            trials_file.close();
        }
    }

    void solveTestTasks(const vector<typename OptMethodType::Task> &tasks, const string &saveDirectory) {
        int numberTasks = (int)tasks.size();

        double totalStartTime = omp_get_wtime();
        for (int i = 0; i < numberTasks; i++) {
            solveTestTask(tasks[i], saveDirectory);
            std::cout << "\n\n";
        }
        double totalEndTime = omp_get_wtime();

        std::cout << "Total time: " << totalEndTime - totalStartTime << "\n";
    }
};

#endif // SOLVER_H_
