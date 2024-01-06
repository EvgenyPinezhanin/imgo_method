#include <iostream>
#include <vector>
#include <string>

#include <Solver.h>
#include <opt_methods/GsaMethod.h>
#include <opt_problems/OneDimensionalProblem.h>
#include <test_opt_problems/OneDimensionalTestProblems.h>
#include <gnuplot/Script.h>

#define CALC
#define DRAW

using Task = opt::Task<OneDimensionalProblem>;
using Parameters = GsaMethod<OneDimensionalProblem>::Parameters;

const std::string methodName = "gsa";
const int displayType = 0; // 0 - application, 1 - png, 2 - png(notitle)
const int problemBlock = 1; // 0 - sample, 1 - test
const int problemNumber = 5; // 0, 1, 2, ...

int main() {
    double reliability = 2.0, accuracy = 0.001;
    int maxTrials = 100000, maxFevals = 100000;
    Parameters parameters(accuracy, 0.0, maxTrials, maxFevals, reliability);

    std::vector<Task> tasks = { Task("Sample Task №1", sampleTasks[0], parameters),
                                Task("Sample Task №2", sampleTasks[1], parameters),
                                Task("Sample Task №3", sampleTasks[2], parameters),
                                Task("Sample Task №4", sampleTasks[3], parameters),

                                Task(  "Test Task №1",   testTasks[0], parameters),
                                Task(  "Test Task №2",   testTasks[1], parameters),
                                Task(  "Test Task №3",   testTasks[2], parameters),
                                Task(  "Test Task №4",   testTasks[3], parameters),
                                Task(  "Test Task №5",   testTasks[4], parameters),
                                Task(  "Test Task №6",   testTasks[5], parameters),
                                Task(  "Test Task №7",   testTasks[6], parameters),
                                Task(  "Test Task №8",   testTasks[7], parameters),
                                Task(  "Test Task №9",   testTasks[8], parameters),
                                Task( "Test Task №10",   testTasks[9], parameters),
                                Task( "Test Task №11",  testTasks[10], parameters),
                                Task( "Test Task №12",  testTasks[11], parameters),
                                Task( "Test Task №13",  testTasks[12], parameters),
                                Task( "Test Task №14",  testTasks[13], parameters),
                                Task( "Test Task №15",  testTasks[14], parameters),
                                Task( "Test Task №16",  testTasks[15], parameters),
                                Task( "Test Task №17",  testTasks[16], parameters),
                                Task( "Test Task №18",  testTasks[17], parameters),
                                Task( "Test Task №19",  testTasks[18], parameters),
                                Task( "Test Task №20",  testTasks[19], parameters) };

#if defined( CALC )
    GsaMethod<OneDimensionalProblem> method;
    Solver solver(method);
    solver.solveTasks(tasks, "output_data/" + methodName + "_test/");
#endif

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
