#include <iostream>
#include <vector>
#include <string>

#include <Solver.h>
#include <opt_methods/ImgoMethod.h>
#include <opt_problems/OneDimConstrainedProblem.h>
#include <test_opt_problems/OneDimConstrainedTestProblems.h>
#include <gnuplot/Script.h>

#define CALC
// #define DRAW

using Task = ImgoMethod<OneDimConstrainedProblem>::Task;
using Parameters = ImgoMethod<OneDimConstrainedProblem>::Parameters;

const std::string methodName = "imgo";
const int displayType = 1; // 0 - application, 1 - png, 2 - png(notitle)
const int functionBlock = 0; // 0 - sample, 1 - test
const int functionNumber = 0; // 0 - f1, 1 - f2, ...

int main() {
    double reliability = 2.0, accuracy = 0.0001, d = 0.0;
    int maxTrials = 100000, maxFevals = 100000;
    Parameters parameters(accuracy, 0.0, maxTrials, maxFevals, reliability, d);

    vector<Task> tasks = { Task("Sample Task №1", blockNames[0],  1, sampleTasks[0], parameters),
                           Task("Sample Task №2", blockNames[0],  2, sampleTasks[1], parameters),
                           Task("Sample Task №3", blockNames[0],  3, sampleTasks[2], parameters),

                           Task(  "Test Task №1", blockNames[1],  1,   testTasks[0], parameters),
                           Task(  "Test Task №2", blockNames[1],  2,   testTasks[1], parameters),
                           Task(  "Test Task №3", blockNames[1],  3,   testTasks[2], parameters),
                           Task(  "Test Task №4", blockNames[1],  4,   testTasks[3], parameters),
                           Task(  "Test Task №5", blockNames[1],  5,   testTasks[4], parameters),
                           Task(  "Test Task №6", blockNames[1],  6,   testTasks[5], parameters),
                           Task(  "Test Task №7", blockNames[1],  7,   testTasks[6], parameters),
                           Task(  "Test Task №8", blockNames[1],  8,   testTasks[7], parameters),
                           Task(  "Test Task №9", blockNames[1],  9,   testTasks[8], parameters),
                           Task( "Test Task №10", blockNames[1], 10,   testTasks[9], parameters) };

#if defined( CALC )
    ImgoMethod<OneDimConstrainedProblem> method;
    Solver solver(method);
    solver.solveTasks(tasks, "output_data/" + methodName + "_test/");
#endif


    int taskArraySize = taskArray.size();
    for (int i = 0; i < taskArraySize; i++) {
        if (taskArray[i].used) {
            imgo.setF(taskArray[i].f);
            imgo.setNumberConstraints(taskArray[i].numberConstraints);
            imgo.setAB(taskArray[i].A[0], taskArray[i].B[0]);
            imgo.setEps(taskArray[i].eps);
            imgo.setMaxTrials(taskArray[i].maxTrials);
            imgo.setMaxFevals(taskArray[i].maxFevals);
            imgo.setR(taskArray[i].r);
            imgo.setD(taskArray[i].d);

            imgo.solve(numberTrials, numberFevals, x);
            imgo.getL(L);

            printResultImgo(taskArray[i].name, taskArray[i].numberConstraints, taskArray[i].A[0], taskArray[i].B[0], taskArray[i].L,
                            taskArray[i].XOpt[0], taskArray[i].f(taskArray[i].XOpt[0], taskArray[i].numberConstraints + 1),
                            taskArray[i].maxTrials, taskArray[i].maxFevals, taskArray[i].eps, taskArray[i].r, taskArray[i].d,
                            numberTrials, numberFevals, L, x, taskArray[i].f(x, taskArray[i].numberConstraints + 1));

            addPointGnuplot(ofstr, taskArray[i].XOpt[0], taskArray[i].f(taskArray[i].XOpt[0], taskArray[i].numberConstraints + 1));
            addPointGnuplot(ofstr, x, taskArray[i].f(x, taskArray[i].numberConstraints + 1));

            imgo.getTrialPoints(trials);
            addPointsGnuplot(ofstr, trials);
        }
    }
    ofstr.close();

    initArrayGnuplot(ofstrOpt, "A", taskArraySize);
    initArrayGnuplot(ofstrOpt, "B", taskArraySize);
    initArrayGnuplot(ofstrOpt, "numberConstraints", taskArraySize);
    for (int i = 0; i < taskArraySize; i++) {
        setValueInArrayGnuplot(ofstrOpt, "A", i + 1, taskArray[i].A[0], false);
        setValueInArrayGnuplot(ofstrOpt, "B", i + 1, taskArray[i].B[0], false);
        setValueInArrayGnuplot(ofstrOpt, "numberConstraints", i + 1, taskArray[i].numberConstraints, false);
    }
    initArrayGnuplot(ofstrOpt, "functionBlockName", taskArraySize);
    initArrayGnuplot(ofstrOpt, "functionNumber", taskArraySize);
    int index = 1;
    for (int i = 0; i < numberBlocks; i++) {
        for (int j = 0; j < numberFunctions[i]; j++) {
            setValueInArrayGnuplot(ofstrOpt, "functionNumber", index, j + 1, false);
            setValueInArrayGnuplot(ofstrOpt, "functionBlockName", index, functionBlockName[i]);
            index++;
        }
    }
    ofstrOpt.close();

    int totalFunctionNumber = functionNumber;
    for (int i = 0; i < functionBlock; i++)
        totalFunctionNumber += numberFunctions[i];
    vector<int> args{ displayType, totalFunctionNumber };
    drawGraphGnuplot("scripts/imgo_test.gp", args);

#if defined( DRAW )
    Script script("scripts/onedimensional_test.gp");
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
