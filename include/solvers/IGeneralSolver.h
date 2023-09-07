#ifndef I_GENERAL_SOLVER_H_
#define I_GENERAL_SOLVER_H_

#include <iostream>
#include <vector>
#include <string>

#include <gnuplot/TrialsFile.h>
#include <omp.h>

template <typename OptMethodType, typename ParametersMethod, typename TaskType,
          typename ResultMethodType, typename ReportType, typename TrialType, typename PointType>
class IGeneralSolver {
protected:
    OptMethodType optMethod;
    ReportType report;

public:
    IGeneralSolver():
        optMethod(),
        report()
    {};

    void solve(const vector<TaskType> &tasks, const string &saveDirectory) {
        TrialsFile trialsFile;
        ReportType report;
        ResultMethodType result;
        vector<TrialType> trials;
        vector<PointType> optimalPoints;
        PointType point;
        int numberTasks = (int)tasks.size();
        double startTime, endTime, workTime;

        double totalStartTime = omp_get_wtime();
        for (int i = 0; i < numberTasks; i++) {
            if (tasks[i].use) {
                optMethod.setProblem(tasks[i].optProblem);
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

                tasks[i].optProblem.getOptimalPoints(optimalPoints);
                trialsFile.addPoints(optimalPoints, tasks[i].optProblem.getOptimalValue());

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

#endif // I_GENERAL_SOLVER_H_
