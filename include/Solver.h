#ifndef SOLVER_H_
#define SOLVER_H_

#include <iostream>
#include <vector>
#include <string>

#include <test_opt_problems/onedimensional_opt_problems.h>
#include <gnuplot/TrialsFile.h>
#include <omp.h>

template <typename OptMethodType>
class Solver {
protected:
    OptMethodType optMethod;
    typename OptMethodType::Parameters parameters;
    typename OptMethodType::Result result;
    typename OptMethodType::Report report;

public:
    Solver(const OptMethodType &_optMethod):
        optMethod(_optMethod),
        result(),
        report()
    {};

    void solve(const vector<typename OptMethodType::Task> &tasks, const string &saveDirectory) {
        TrialsFile trialsFile;

        vector<typename OptMethodType::Trial> trials;
        vector<typename OptMethodType::OptProblem::Point> optimalPoints;
        typename OptMethodType::OptProblem::Point point;
        int numberTasks = (int)tasks.size();
        double startTime, endTime, workTime;

        double totalStartTime = omp_get_wtime();
        for (int i = 0; i < numberTasks; i++) {
            if (tasks[i].use) {
                optMethod.setProblem(tasks[i].problem);
                optMethod.setParameters(tasks[i].parameters);

                startTime = omp_get_wtime();
                optMethod.solve(result);
                endTime = omp_get_wtime();
                workTime = endTime - startTime;

                report.print(std::cout, tasks[i], result, workTime);
                std::cout << "\n\n";

                trialsFile.newFile(saveDirectory + blockNames[tasks[i].blockNumber] + 
                                   "_" + std::to_string(tasks[i].functionNumber));
                if (!trialsFile.isOpen()) std::cerr << "TrialsFile opening error\n";

                tasks[i].problem.getOptimalPoints(optimalPoints);
                trialsFile.addPoints(optimalPoints, tasks[i].problem.getOptimalValue());

                trialsFile.addPoint(result.point, result.value);

                optMethod.getTrialPoints(trials);
                trialsFile.addPoints(trials);
            }
        }
        double totalEndTime = omp_get_wtime();
        std::cout << "Total time: " << totalEndTime - totalStartTime << "\n";

        trialsFile.closeFile();
    }
};

#endif // SOLVER_H_
