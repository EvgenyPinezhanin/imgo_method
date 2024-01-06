#include <iostream>
#include <vector>
#include <string>

#include <Solver.h>
#include <opt_methods/ImgoMethod.h>
#include <opt_problems/OneDimensionalConstrainedProblem.h>
#include <test_opt_problems/OneDimensionalConstrainedTestProblems.h>
#include <gnuplot/Script.h>

#define CALC
#define DRAW

using Task = opt::Task<OneDimensionalConstrainedProblem>;
using Parameters = ImgoMethod<OneDimensionalConstrainedProblem>::Parameters;

const std::string methodName = "imgo";
const int displayType = 1; // 0 - application, 1 - png, 2 - png(notitle)
const int problemBlock = 0; // 0 - sample, 1 - test
const int problemNumber = 2; // 0, 1, 2, ...

int main() {
    double reliability = 2.0, accuracy = 0.0001, d = 0.0;
    int maxTrials = 100000, maxFevals = 100000;
    std::vector<Parameters> parameters = {
        Parameters(accuracy, 0.0, maxTrials, maxFevals, std::vector<double>{ reliability }, d),
        Parameters(accuracy, 0.0, maxTrials, maxFevals, std::vector<double>{ reliability, reliability }, d),
        Parameters(accuracy, 0.0, maxTrials, maxFevals, std::vector<double>{ reliability, reliability, reliability }, d),

        Parameters(accuracy, 0.0, maxTrials, maxFevals, std::vector<double>{ reliability, reliability }, d),
        Parameters(accuracy, 0.0, maxTrials, maxFevals, std::vector<double>{ reliability, reliability }, d),
        Parameters(accuracy, 0.0, maxTrials, maxFevals, std::vector<double>{ reliability, reliability }, d),
        Parameters(accuracy, 0.0, maxTrials, maxFevals, std::vector<double>{ reliability, reliability, reliability }, d),
        Parameters(accuracy, 0.0, maxTrials, maxFevals, std::vector<double>{ reliability, reliability, reliability }, d),
        Parameters(accuracy, 0.0, maxTrials, maxFevals, std::vector<double>{ reliability, reliability, reliability }, d),
        Parameters(accuracy, 0.0, maxTrials, maxFevals, std::vector<double>{ reliability, reliability, reliability }, d),
        Parameters(accuracy, 0.0, maxTrials, maxFevals, std::vector<double>{ reliability, reliability, reliability, reliability }, d),
        Parameters(accuracy, 0.0, maxTrials, maxFevals, std::vector<double>{ reliability, reliability, reliability, reliability }, d),
        Parameters(accuracy, 0.0, maxTrials, maxFevals, std::vector<double>{ reliability, reliability, reliability, reliability }, d)
    };

    std::vector<Task> tasks = { Task("Sample Task №1", sampleTasks[0],  parameters[0]),
                                Task("Sample Task №2", sampleTasks[1],  parameters[1]),
                                Task("Sample Task №3", sampleTasks[2],  parameters[2]),

                                Task(  "Test Task №1", testTasks[0],  parameters[3]),
                                Task(  "Test Task №2", testTasks[1],  parameters[4]),
                                Task(  "Test Task №3", testTasks[2],  parameters[5]),
                                Task(  "Test Task №4", testTasks[3],  parameters[6]),
                                Task(  "Test Task №5", testTasks[4],  parameters[7]),
                                Task(  "Test Task №6", testTasks[5],  parameters[8]),
                                Task(  "Test Task №7", testTasks[6],  parameters[9]),
                                Task(  "Test Task №8", testTasks[7], parameters[10]),
                                Task(  "Test Task №9", testTasks[8], parameters[11]),
                                Task( "Test Task №10", testTasks[9], parameters[12]) };

#if defined( CALC )
    ImgoMethod<OneDimensionalConstrainedProblem> method;
    Solver solver(method);
    solver.solveTasks(tasks, "output_data/" + methodName + "_test/");
#endif

#if defined( DRAW )
    Script script("scripts/onedimensional_constrained_test.gp");
    script.addArg(methodName);
    script.addArgs(std::vector<int>{ displayType, problemBlock, problemNumber });
    script.start();
    if (script.isError() == 2) std::cerr << "Error gnuplot\n";
    if (script.isError() == 1) std::cerr << "Error chmod\n";
#endif

#if defined( _MSC_VER )
    cin.get();
#endif

    return 0;
}
