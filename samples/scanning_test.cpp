#include <iostream>
#include <vector>
#include <string>

#include <opt_methods/ScanningMethod.h>
#include <tasks/ScanningTask.h>
#include <test_opt_problems/onedimensional_opt_problems.h>
#include <print_result.h>
#include <gnuplot/TrialsFile.h>
#include <gnuplot/Script.h>
#include <omp.h>

using std::string;
using std::to_string;
using std::vector;
using std::cout;
using std::cerr;

#define CALC

const string methodName = "scanning";
const int displayType = 1; // 0 - application, 1 - png, 2 - png(notitle)
const int functionBlock = 1; // 0 - sample, 1 - test
const int functionNumber = 5; // 0 - f1, 1 - f2, ...

int main() {
    double accuracy = 0.001;
    int maxTrials = 100000, maxFevals = 100000;

    vector<ScanningTask> tasks = { ScanningTask("Sample Task №1", sampleTasks[0], 0, 1,  accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Sample Task №2", sampleTasks[1], 0, 2,  accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Sample Task №3", sampleTasks[2], 0, 3,  accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Sample Task №4", sampleTasks[3], 0, 4,  accuracy, 0.0, maxTrials, maxFevals),

                                   ScanningTask("Test Task №1",   testTasks[0],   1, 1,  accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №2",   testTasks[1],   1, 2,  accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №3",   testTasks[2],   1, 3,  accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №4",   testTasks[3],   1, 4,  accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №5",   testTasks[4],   1, 5,  accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №6",   testTasks[5],   1, 6,  accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №7",   testTasks[6],   1, 7,  accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №8",   testTasks[7],   1, 8,  accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №9",   testTasks[8],   1, 9,  accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №10",  testTasks[9],   1, 10, accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №11",  testTasks[10],  1, 11, accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №12",  testTasks[11],  1, 12, accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №13",  testTasks[12],  1, 13, accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №14",  testTasks[13],  1, 14, accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №15",  testTasks[14],  1, 15, accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №16",  testTasks[15],  1, 16, accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №17",  testTasks[16],  1, 17, accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №18",  testTasks[17],  1, 18, accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №19",  testTasks[18],  1, 19, accuracy, 0.0, maxTrials, maxFevals),
                                   ScanningTask("Test Task №20",  testTasks[19],  1, 20, accuracy, 0.0, maxTrials, maxFevals) };

    ScanningMethod scanning;

    TrialsFile trialsFile;
    ResultMethod result;
    vector<Trial> trials;
    vector<double> optimalPoints;
    double point;
    int numberTasks = (int)tasks.size();
    double startTime, endTime, workTime;

    double totalStartTime = omp_get_wtime();
#if defined( CALC )
    for (int i = 0; i < numberTasks; i++) {
        if (tasks[i].use) {
            scanning.setProblem(tasks[i].optProblem);
            scanning.setAccuracy(tasks[i].accuracy);
            scanning.setMaxTrials(tasks[i].maxTrials);
            scanning.setMaxFevals(tasks[i].maxFevals);

            startTime = omp_get_wtime();
            scanning.solve(result);
            endTime = omp_get_wtime();
            workTime = endTime - startTime;

            printResult(cout, tasks[i], result, workTime);

            trialsFile.newFile("output_data/" + methodName + "_test/" + blockNames[tasks[i].blockNumber] + 
                               "_" + to_string(tasks[i].functionNumber));
            if (!trialsFile.isOpen()) cerr << "TrialsFile opening error\n";

            tasks[i].optProblem.getOptimalPoints(optimalPoints);
            trialsFile.addPoints(optimalPoints, tasks[i].optProblem.computeObjFunction(optimalPoints[0]));
            point = result.solution;
            trialsFile.addPoint(point, tasks[i].optProblem.computeObjFunction(point));

            scanning.getTrialPoints(trials);
            trialsFile.addPoints(trials);
        }
    }
    trialsFile.closeFile();
#endif
    double totalEndTime = omp_get_wtime();
    cout << "Total time: " << totalEndTime - totalStartTime << "\n";

    Script script("scripts/onedimensional_test.gp");
    script.addArg(methodName);
    script.addArgs(vector<int>{ displayType, functionBlock, functionNumber });
    script.start();
    if (script.isError() == 1) cerr << "Error chmod\n";
    if (script.isError() == 2) cerr << "Error gnuplot\n";

#if defined( _MSC_VER )
    cin.get();
#endif

    return 0;
}
