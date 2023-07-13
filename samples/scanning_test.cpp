#include <iostream>
#include <vector>

#include <opt_methods/ScanningMethod.h>
#include <tasks/ScanningTask.h>
#include <test_opt_problems/one_dimensional/sample_opt_problems.h>
#include <test_opt_problems/one_dimensional/test_opt_problems.h>
#include <print_result.h>
#include <gnuplot/TrialsFile.h>
#include <gnuplot/VariablesFile.h>
#include <gnuplot/start_script.h>

using namespace std;

const int functionNumber = 0; // 0 - f1, 1 - f2, ...
const int functionBlock = 0; // 0 - sample, 1 - test
const int displayType = 1; // 0 - application, 1 - png, 2 - png(notitle)

const int numberBlocks = 2;
const int numberFunctions[2] = { 4, 20 };
const vector<string> functionBlockName{ "sample", "test" };

int main() {
    TrialsFile trialsFile("output_data/scanning_test.txt");
    if (!trialsFile.isOpen()) cerr << "TrialsFile opening error\n";
    VariablesFile variablesFile("output_data/scanning_test_vars.txt");
    if (!variablesFile.isOpen()) cerr << "File opening error\n";

    double accuracy = 0.001;
    int maxTrials = 100000, maxFevals = 100000;

    vector<ScanningTask> tasks = { ScanningTask("Sample Task №1", sampleTasks[0], accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Sample Task №2", sampleTasks[1], accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Sample Task №3", sampleTasks[2], accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Sample Task №4", sampleTasks[3], accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №1", testTasks[0], accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №2", testTasks[1], accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №3", testTasks[2], accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №4", testTasks[3], accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №5", testTasks[4], accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №6", testTasks[5], accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №7", testTasks[6], accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №8", testTasks[7], accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №9", testTasks[8], accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №10", testTasks[9], accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №11", testTasks[10], accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №12", testTasks[11], accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №13", testTasks[12], accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №14", testTasks[13], accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №15", testTasks[14], accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №16", testTasks[15], accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №17", testTasks[16], accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №18", testTasks[17], accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №19", testTasks[18], accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №20", testTasks[19], accuracy, 0.0, maxTrials, maxFevals) };

    ScanningMethod scanning;

    ResultMethod result;
    vector<Trial> trials;
    double optimalPoint, point;

    int numberTasks = (int)tasks.size();
    for (int i = 0; i < numberTasks; i++) {
        if (tasks[i].use) {
            scanning.setProblem(tasks[i].optProblem);
            scanning.setAccuracy(tasks[i].accuracy);
            scanning.setMaxTrials(tasks[i].maxTrials);
            scanning.setMaxFevals(tasks[i].maxFevals);

            scanning.solve(result);

            printResult(cout, tasks[i], result);

            optimalPoint = tasks[i].optProblem.getOptimalPoint();
            trialsFile.addPoint(optimalPoint, tasks[i].optProblem.computeObjFunction(optimalPoint));
            point = result.solution;
            trialsFile.addPoint(point, tasks[i].optProblem.computeObjFunction(point));

            scanning.getTrialPoints(trials);
            trialsFile.addPoints(trials);
        }
    }
    trialsFile.closeFile();

    variablesFile.initArray("A", numberTasks);
    variablesFile.initArray("B", numberTasks);
    for (int i = 0; i < numberTasks; i++) {
        variablesFile.setValueInArray("A", i + 1, tasks[i].optProblem.getSearchArea().getLowerBound(), false);
        variablesFile.setValueInArray("B", i + 1, tasks[i].optProblem.getSearchArea().getUpBound(), false);
    }
    variablesFile.initArray("functionBlockName", numberTasks);
    variablesFile.initArray("functionNumber", numberTasks);
    int index = 1;
    for (int i = 0; i < numberBlocks; i++) {
        for (int j = 0; j < numberFunctions[i]; j++) {
            variablesFile.setValueInArray("functionNumber", index, j + 1, false);
            variablesFile.setValueInArray("functionBlockName", index, functionBlockName[i]);
            index++;
        }
    }
    variablesFile.closeFile();

    int totalFunctionNumber = functionNumber;
    for (int i = 0; i < functionBlock; i++)
        totalFunctionNumber += numberFunctions[i];
    vector<int> args{ displayType, totalFunctionNumber };
    startScript("scripts/gsa_test.gp", args);
    if (errorScript == 1) cerr << "Error chmod" << endl;
    if (errorScript == 2) cerr << "Error gnuplot" << endl;

#if defined( _MSC_VER )
    cin.get();
#endif

    return 0;
}
