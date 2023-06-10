#if defined( _MSC_VER )
    #define _CRT_SECURE_NO_WARNINGS    
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <limits>

#include <mggsa.h>
#include <task.h>
#include <output_results.h>

using namespace std;

const int functionNumber = 2; // 1 - f1, 2 - f2, 3 - f3, 4 - f4
const int displayType = 1; // 0 - application, 1 - png, 2 - png(notitle)

double f1(vector<double> x, int j) {
    switch (j) {
        case 1: return 1.0 - x[0] - x[1];
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f2(vector<double> x, int j) {
    switch (j) {
        case 1: return (x[0] - 1.0) * (x[0] - 1.0) / 5.0 + (x[1] - 1.0) * (x[1] - 1.0) / 5.0;
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f3(vector<double> x, int j) {
    switch (j) {
        case 1: return 1.0 - x[0] - x[1];
        case 2: return x[0] * x[0] / 5.0 + x[1] * x[1] / 5.0;
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f4(vector<double> x, int j) {
    switch (j) {
        case 1: return (x[0] - 2.0) * (x[0] - 2.0) + (x[1] - 2.0) * (x[1] - 2.0) - 2.0;
        case 2: return x[0] * x[0] / 5.0 + x[1] * x[1] / 5.0;
        default: return numeric_limits<double>::quiet_NaN();
    }
}

int main() {
    ofstream ofstr("output_data/mggsa_sample.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstrOpt("output_data/mggsa_sample_opt.txt");
    if (!ofstrOpt.is_open()) cerr << "File opening error\n";

    double eps = 0.01, r = 2.2, d = 0.0;
    int n = 2, den = 10, key = 1;
    int maxTrials = 100000, maxFevals = 100000;

    vector<TaskMggsa> taskArray = { TaskMggsa(f1, "f1", n, 0, vector<double>{-4.0, -4.0}, vector<double>{4.0, 4.0},
                                              vector<double>{4.0, 4.0}, vector<double>{}, eps, maxTrials, maxFevals, r, d, den, key, -1),
                                    TaskMggsa(f2, "f2", n, 0, vector<double>{-4.0, -4.0}, vector<double>{4.0, 4.0},
                                              vector<double>{1.0, 1.0}, vector<double>{}, eps, maxTrials, maxFevals, r, d, den, key, -1),
                                    TaskMggsa(f3, "f3", n, 1, vector<double>{-1.0, -1.0}, vector<double>{1.0, 1.0},
                                              vector<double>{0.5, 0.5}, vector<double>{}, eps, maxTrials, maxFevals, r, d, den, key, -1),
                                    TaskMggsa(f4, "f4", n, 1, vector<double>{0.0, 0.0}, vector<double>{3.0, 3.0},
                                              vector<double>{1.0, 1.0}, vector<double>{}, eps, maxTrials, maxFevals, r, d, den, key, -1) };

    MggsaMethod mggsa;

    vector<double> X, L;
    vector<vector<double>> points;
    vector<TrialConstrained> trials;
    int numberTrials, numberFevals;

    for (int i = 0; i < taskArray.size(); i++) {
        if (taskArray[i].used) {
            mggsa.setF(taskArray[i].f);
            mggsa.setN(taskArray[i].n);
            mggsa.setNumberConstraints(taskArray[i].numberConstraints);
            mggsa.setAB(taskArray[i].A, taskArray[i].B);
            mggsa.setEps(taskArray[i].eps);
            mggsa.setMaxTrials(taskArray[i].maxTrials);
            mggsa.setMaxFevals(taskArray[i].maxFevals);
            mggsa.setR(taskArray[i].r);
            mggsa.setD(taskArray[i].d);
            mggsa.setDen(taskArray[i].den);
            mggsa.setKey(taskArray[i].key);

            mggsa.solve(numberTrials, numberFevals, X);
            mggsa.getL(L);

            printResultMggsa(taskArray[i].name, taskArray[i].n, taskArray[i].numberConstraints, taskArray[i].A, taskArray[i].B,
                             taskArray[i].L, taskArray[i].XOpt, taskArray[i].f(taskArray[i].XOpt, taskArray[i].numberConstraints + 1),
                             taskArray[i].maxTrials, taskArray[i].maxFevals, taskArray[i].eps, taskArray[i].r, taskArray[i].d,
                             taskArray[i].den, taskArray[i].key, taskArray[i].incr, numberTrials, numberFevals, L, X,
                             taskArray[i].f(X, taskArray[i].numberConstraints + 1));

            addPointGnuplot(ofstr, taskArray[i].XOpt, taskArray[i].f(taskArray[i].XOpt, taskArray[i].numberConstraints + 1));
            addPointGnuplot(ofstr, X, taskArray[i].f(X, taskArray[i].numberConstraints + 1));

            mggsa.getPoints(points);
            mggsa.getTrialPoints(trials);
            addPointsGnuplot(ofstr, points, trials);
        }
    }
    ofstr.close();

    size_t sizeTaskArray = taskArray.size();
    string arraysName[] = { "AX", "BX", "AY", "BY" };
    for (int i = 0; i < sizeTaskArray; i++) {
        initArrayGnuplot(ofstrOpt, arraysName[i], sizeTaskArray);
    }
    for (int i = 0; i < sizeTaskArray; i++) {
        for (int j = 0; j < 2; j++) {
            setValueInArrayGnuplot(ofstrOpt, arraysName[2 * j], i + 1, taskArray[i].A[j], false);
            setValueInArrayGnuplot(ofstrOpt, arraysName[2 * j + 1], i + 1, taskArray[i].B[j], false);
        }
    }
    ofstrOpt.close();

    vector<int> args{ displayType, functionNumber };
    drawGraphGnuplot("scripts/mggsa_sample.gp", args);

#if defined(_MSC_VER)
    cin.get();
#endif

    return 0;
}
