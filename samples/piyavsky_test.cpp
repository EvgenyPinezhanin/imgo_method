#include <iostream>
#include <vector>
#include <string>

#include <opt_methods/PiyavskyMethod.h>
#include <tasks/PiyavskyGsaTask.h>
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

const string methodName = "piyavsky";
const int displayType = 1; // 0 - application, 1 - png, 2 - png(notitle)
const int functionBlock = 1; // 0 - sample, 1 - test
const int functionNumber = 5; // 0 - f1, 1 - f2, ...

int main() {
    double reliability = 2.0, accuracy = 0.001;
    int maxTrials = 100000, maxFevals = 100000;

    vector<PiyavskyGsaTask> tasks = { PiyavskyGsaTask("Sample Task №1", sampleTasks[0], 0, 1,  reliability, accuracy,
                                                      0.0, maxTrials, maxFevals),
                                      PiyavskyGsaTask("Sample Task №2", sampleTasks[1], 0, 2,  reliability, accuracy,
                                                      0.0, maxTrials, maxFevals),
                                      PiyavskyGsaTask("Sample Task №3", sampleTasks[2], 0, 3,  reliability, accuracy,
                                                      0.0, maxTrials, maxFevals),
                                      PiyavskyGsaTask("Sample Task №4", sampleTasks[3], 0, 4,  reliability, accuracy,
                                                      0.0, maxTrials, maxFevals),

                                      PiyavskyGsaTask("Test Task №1",   testTasks[0],   1, 1,  reliability, accuracy,
                                                      0.0, maxTrials, maxFevals),
                                      PiyavskyGsaTask("Test Task №2",   testTasks[1],   1, 2,  reliability, accuracy,
                                                      0.0, maxTrials, maxFevals),
                                      PiyavskyGsaTask("Test Task №3",   testTasks[2],   1, 3,  reliability, accuracy,
                                                      0.0, maxTrials, maxFevals),
                                      PiyavskyGsaTask("Test Task №4",   testTasks[3],   1, 4,  reliability, accuracy,
                                                      0.0, maxTrials, maxFevals),
                                      PiyavskyGsaTask("Test Task №5",   testTasks[4],   1, 5,  reliability, accuracy,
                                                      0.0, maxTrials, maxFevals),
                                      PiyavskyGsaTask("Test Task №6",   testTasks[5],   1, 6,  reliability, accuracy,
                                                      0.0, maxTrials, maxFevals),
                                      PiyavskyGsaTask("Test Task №7",   testTasks[6],   1, 7,  reliability, accuracy,
                                                      0.0, maxTrials, maxFevals),
                                      PiyavskyGsaTask("Test Task №8",   testTasks[7],   1, 8,  reliability, accuracy,
                                                      0.0, maxTrials, maxFevals),
                                      PiyavskyGsaTask("Test Task №9",   testTasks[8],   1, 9,  reliability, accuracy,
                                                      0.0, maxTrials, maxFevals),
                                      PiyavskyGsaTask("Test Task №10",  testTasks[9],   1, 10, reliability, accuracy,
                                                      0.0, maxTrials, maxFevals),
                                      PiyavskyGsaTask("Test Task №11",  testTasks[10],  1, 11, reliability, accuracy,
                                                      0.0, maxTrials, maxFevals),
                                      PiyavskyGsaTask("Test Task №12",  testTasks[11],  1, 12, reliability, accuracy,
                                                      0.0, maxTrials, maxFevals),
                                      PiyavskyGsaTask("Test Task №13",  testTasks[12],  1, 13, reliability, accuracy,
                                                      0.0, maxTrials, maxFevals),
                                      PiyavskyGsaTask("Test Task №14",  testTasks[13],  1, 14, reliability, accuracy,
                                                      0.0, maxTrials, maxFevals),
                                      PiyavskyGsaTask("Test Task №15",  testTasks[14],  1, 15, reliability, accuracy,
                                                      0.0, maxTrials, maxFevals),
                                      PiyavskyGsaTask("Test Task №16",  testTasks[15],  1, 16, reliability, accuracy,
                                                      0.0, maxTrials, maxFevals),
                                      PiyavskyGsaTask("Test Task №17",  testTasks[16],  1, 17, reliability, accuracy,
                                                      0.0, maxTrials, maxFevals),
                                      PiyavskyGsaTask("Test Task №18",  testTasks[17],  1, 18, reliability, accuracy,
                                                      0.0, maxTrials, maxFevals),
                                      PiyavskyGsaTask("Test Task №19",  testTasks[18],  1, 19, reliability, accuracy,
                                                      0.0, maxTrials, maxFevals),
                                      PiyavskyGsaTask("Test Task №20",  testTasks[19],  1, 20, reliability, accuracy,
                                                      0.0, maxTrials, maxFevals) };

    PiyavskyMethod piyavsky;

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
            piyavsky.setProblem(tasks[i].optProblem);
            piyavsky.setReliability(tasks[i].reliability);
            piyavsky.setAccuracy(tasks[i].accuracy);
            piyavsky.setMaxTrials(tasks[i].maxTrials);
            piyavsky.setMaxFevals(tasks[i].maxFevals);

            startTime = omp_get_wtime();
            piyavsky.solve(result);
            endTime = omp_get_wtime();
            workTime = endTime - startTime;

            printResult(cout, tasks[i], result, workTime);

            trialsFile.newFile("output_data/" + methodName + "_test/" + blockNames[tasks[i].blockNumber] + 
                               "_" + to_string(tasks[i].functionNumber));
            if (!trialsFile.isOpen()) cerr << "TrialsFile opening error\n";

            tasks[i].optProblem.getOptimalPoints(optimalPoints);
            trialsFile.addPoints(optimalPoints, tasks[i].optProblem.getOptimalValue());

            trialsFile.addPoint(result.point, result.value);

            piyavsky.getTrialPoints(trials);
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
