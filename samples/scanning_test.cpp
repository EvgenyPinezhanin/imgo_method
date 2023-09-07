#include <iostream>
#include <string>
#include <vector>

#include <tasks/ScanningTask.h>
#include <test_opt_problems/onedimensional_opt_problems.h>
#include <solvers/ScanningSolver.h>
#include <gnuplot/Script.h>

using std::cerr;
using std::string;
using std::vector;

#define CALC

const string methodName = "scanning";
const bool isDraw = false;
const int displayType = 1; // 0 - application, 1 - png, 2 - png(notitle)
const int problemBlock = 1; // 0 - sample, 1 - test
const int problemNumber = 5; // 0, 1, 2, ...

int main() {
    double accuracy = 0.001;
    int maxTrials = 100000, maxFevals = 100000;

    vector<ScanningTask> tasks = { ScanningTask("Sample Task №1", 0,  1, sampleTasks[0],
                                                ScanningParameters(accuracy, 0.0, maxTrials, maxFevals)),
                                   ScanningTask("Sample Task №2", 0,  2, sampleTasks[1],
                                                ScanningParameters(accuracy, 0.0, maxTrials, maxFevals)),
                                   ScanningTask("Sample Task №3", 0,  3, sampleTasks[2],
                                                ScanningParameters(accuracy, 0.0, maxTrials, maxFevals)),
                                   ScanningTask("Sample Task №4", 0,  4, sampleTasks[3],
                                                ScanningParameters(accuracy, 0.0, maxTrials, maxFevals)),

                                   ScanningTask(  "Test Task №1", 1,  1,   testTasks[0],
                                                ScanningParameters(accuracy, 0.0, maxTrials, maxFevals)),
                                   ScanningTask(  "Test Task №2", 1,  2,   testTasks[1],
                                                ScanningParameters(accuracy, 0.0, maxTrials, maxFevals)),
                                   ScanningTask(  "Test Task №3", 1,  3,   testTasks[2],
                                                ScanningParameters(accuracy, 0.0, maxTrials, maxFevals)),
                                   ScanningTask(  "Test Task №4", 1,  4,   testTasks[3],
                                                ScanningParameters(accuracy, 0.0, maxTrials, maxFevals)),
                                   ScanningTask(  "Test Task №5", 1,  5,   testTasks[4],
                                                ScanningParameters(accuracy, 0.0, maxTrials, maxFevals)),
                                   ScanningTask(  "Test Task №6", 1,  6,   testTasks[5],
                                                ScanningParameters(accuracy, 0.0, maxTrials, maxFevals)),
                                   ScanningTask(  "Test Task №7", 1,  7,   testTasks[6],
                                                ScanningParameters(accuracy, 0.0, maxTrials, maxFevals)),
                                   ScanningTask(  "Test Task №8", 1,  8,   testTasks[7],
                                                ScanningParameters(accuracy, 0.0, maxTrials, maxFevals)),
                                   ScanningTask(  "Test Task №9", 1,  9,   testTasks[8],
                                                ScanningParameters(accuracy, 0.0, maxTrials, maxFevals)),
                                   ScanningTask( "Test Task №10", 1, 10,   testTasks[9],
                                                ScanningParameters(accuracy, 0.0, maxTrials, maxFevals)),
                                   ScanningTask( "Test Task №11", 1, 11,  testTasks[10],
                                                ScanningParameters(accuracy, 0.0, maxTrials, maxFevals)),
                                   ScanningTask( "Test Task №12", 1, 12,  testTasks[11],
                                                ScanningParameters(accuracy, 0.0, maxTrials, maxFevals)),
                                   ScanningTask( "Test Task №13", 1, 13,  testTasks[12],
                                                ScanningParameters(accuracy, 0.0, maxTrials, maxFevals)),
                                   ScanningTask( "Test Task №14", 1, 14,  testTasks[13],
                                                ScanningParameters(accuracy, 0.0, maxTrials, maxFevals)),
                                   ScanningTask( "Test Task №15", 1, 15,  testTasks[14],
                                                ScanningParameters(accuracy, 0.0, maxTrials, maxFevals)),
                                   ScanningTask( "Test Task №16", 1, 16,  testTasks[15],
                                                ScanningParameters(accuracy, 0.0, maxTrials, maxFevals)),
                                   ScanningTask( "Test Task №17", 1, 17,  testTasks[16],
                                                ScanningParameters(accuracy, 0.0, maxTrials, maxFevals)),
                                   ScanningTask( "Test Task №18", 1, 18,  testTasks[17],
                                                ScanningParameters(accuracy, 0.0, maxTrials, maxFevals)),
                                   ScanningTask( "Test Task №19", 1, 19,  testTasks[18],
                                                ScanningParameters(accuracy, 0.0, maxTrials, maxFevals)),
                                   ScanningTask( "Test Task №20", 1, 20,  testTasks[19],
                                                ScanningParameters(accuracy, 0.0, maxTrials, maxFevals)) };

#if defined( CALC )
    ScanningSolver solver;
    solver.solve(tasks, "output_data/" + methodName + "_test/");
#endif

    if (isDraw) {
        Script script("scripts/onedimensional_test.gp");
        script.addArg(methodName);
        script.addArgs(std::vector<int>{ displayType, problemBlock, problemNumber });
        script.start();
        if (script.isError() == 1) std::cerr << "Error chmod\n";
        if (script.isError() == 2) std::cerr << "Error gnuplot\n";
    }

#if defined( _MSC_VER )
    cin.get();
#endif

    return 0;
}
