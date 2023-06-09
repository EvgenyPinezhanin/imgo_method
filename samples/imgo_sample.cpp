#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#if defined( _MSC_VER )
    #define _USE_MATH_DEFINES
    #include <math.h>
#else
    #include <cmath>
#endif

#include <imgo.h>
#include <task.h>
#include <output_results.h>

using namespace std;

const int functionNumber = 1; // 1 - f1, 2 - f2, 3 - f3
const int displayType = 1; // 0 - application, 1 - png

double f1(double x, int j) {
    switch (j) {
        case 1: return sin(x);
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f2(double x, int j) {
    switch(j) {
        case 1: return sin(x);
        case 2: return -2.0 * x + 3.0;
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f3(double x, int j) {
    switch (j) {
        case 1: return x * x - 0.05;
        case 2: return -x + 0.1;
        case 3: return 5.0 * x * x + 3.0 * x - 1.0;
        default: return numeric_limits<double>::quiet_NaN();
    }
}

int main() {
    ofstream ofstr("output_data/imgo_sample.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstrOpt("output_data/imgo_sample_opt.txt");
    if (!ofstrOpt.is_open()) cerr << "File opening error\n";

    double eps = 0.001, r = 3.0, d = 0.0;
    int maxTrials = 100000, maxFevals = 100000;

    vector<TaskImgo> taskArray = { TaskImgo(f1, "f1(x)", 0, -4.0, 4.0, -M_PI / 2.0, vector<double>{ 1.0 }, eps, maxTrials, maxFevals, r, d),
                                   TaskImgo(f2, "f2(x)", 1, 2.0, 8.0, 2.0 * M_PI, vector<double>{ 1.0, 2.0 }, eps, maxTrials, maxFevals, r, d),
                                   TaskImgo(f3, "f3(x)", 2, -2.0, 2.0, 0.1, vector<double>{ 4.0, 1.0, 23.0 }, eps, maxTrials, maxFevals, r, d) };

    ImgoMethod imgo;

    double x;
    int numberTrials, numberFevals;
    vector<double> L;
    vector<TrialConstrained> trials;

    for (int i = 0; i < taskArray.size(); i++) {
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

    initArrayGnuplot(ofstrOpt, "numberConstraints", taskArray.size());
    for (int i = 0; i < taskArray.size(); i++) {
        setValueInArrayGnuplot(ofstrOpt, "numberConstraints", i + 1, taskArray[i].numberConstraints, false);
    }
    ofstrOpt.close();

    vector<int> args{ displayType, functionNumber };
    drawGraphGnuplot("scripts/imgo_sample.gp", args);

#if defined( _MSC_VER )
    cin.get();
#endif

    return 0;
}
