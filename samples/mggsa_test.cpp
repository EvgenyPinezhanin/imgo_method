#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include <mggsa.h>
#include <task.h>
#include <output_results.h>

using namespace std;

const int functionNumber = 1; // 1 - f1, 2 - f2
const int displayType = 2; // 0 - application, 1 - png, 2 - png(notitle)

double f1(vector<double> x, int j) {
    switch (j) {
        case 1: return 0.01 * (pow((x[0] - 2.2), 2) + pow((x[1] - 1.2), 2) - 2.25);
        case 2: return 100.0 * (1.0 - pow((x[0] - 2.0), 2) / 1.44 - pow(0.5 * x[1], 2));
        case 3: return 10.0 * (x[1] - 1.5 - 1.5 * sin(6.283 * (x[0] - 1.75)));
        case 4: return -1.5 * x[0] * x[0] * exp(1.0 - x[0] * x[0] - 20.25 * pow((x[0] - x[1]), 2)) -
                       pow(0.5 * (x[0] - 1.0) * (x[1] - 1.0), 4) * exp(2.0 - pow(0.5 * (x[0] - 1.0), 4) -
                       pow(x[1] - 1.0, 4));
        default: return numeric_limits<double>::quiet_NaN();
    }
}

const double C[20] = { 75.1963666677, -3.8112755343, 0.1269366345, -0.0020567665, 0.000010345,
                      -6.8306567631, 0.0302344793, -0.0012813448, 0.0000352559, -0.0000002266,
                       0.2564581253, -0.0034604030, 0.0000135139, -28.1064434908, -0.0000052375,
                      -0.0000000063, 0.0000000007, 0.0003405462, -0.0000016638, -2.8673112392 };

double f2(vector<double> x, int j) {
    switch (j) {
        case 1: return 450.0 - x[0] * x[1];
        case 2: return (0.1 * x[0] - 1.0) * (0.1 * x[0] - 1.0) - x[1];
        case 3: return 8.0 * (x[0] - 40.0) - (x[1] - 30.0) * (x[1] - 55.0);
        case 4: return x[1] + (x[0] - 35.0) * (x[0] - 30.0) / 125.0 - 80.0;
        case 5: return -(C[0] + C[1] * x[0] + C[2] * x[0] * x[0] + C[3] * pow(x[0], 3) + C[4] * pow(x[0], 4) + C[5] * x[1] +
                         C[6] * x[0] * x[1] + C[7] * x[0] * x[0] * x[1] + C[8] * pow(x[0], 3) * x[1] + C[9] * pow(x[0], 4) * x[1] +
                         C[10] * x[1] * x[1] + C[11] * pow(x[1], 3) + C[12] * pow(x[1], 4) + C[13] / (x[1] + 1) + C[14] * x[0] * x[0] * x[1] * x[1] +
                         C[15] * pow(x[0], 3) * x[1] * x[1] + C[16] * pow(x[0], 3) * pow(x[1], 3) + C[17] * x[0] * x[1] * x[1] +
                         C[18] * x[0] * pow(x[1], 3) + C[19] * exp(0.0005 * x[0] * x[1]));
        default: return numeric_limits<double>::quiet_NaN();
    }
}

int main() {
    ofstream ofstr("output_data/mggsa_test.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstrOpt("output_data/mggsa_test_opt.txt");
    if (!ofstrOpt.is_open()) cerr << "File opening error\n";

    double eps = 0.001, r = 2.2, d = 0.05;
    int n = 2, den = 10, key = 1;
    int maxTrials = 100000, maxFevals = 100000;

    vector<TaskMggsa> taskArray = { TaskMggsa(f1, "f1", n, 3, vector<double>{0.0, -1.0}, vector<double>{4.0, 3.0},
                                              vector<double>{0.942, 0.944}, vector<double>{}, eps, maxTrials, maxFevals,
                                              r, d, 12, key, -1, true),
                                    TaskMggsa(f2, "f2", n, 4, vector<double>{0.0, 0.0}, vector<double>{80.0, 80.0},
                                              vector<double>{77.489, 63.858}, vector<double>{}, eps, maxTrials, maxFevals,
                                              3.3, 0.01, den, key, -1, true) };

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
    vector<double> incrementalParam{  -1.6,  0.3,  1.0,
                                     -50.0, 10.0, 50.0 };
    initArrayGnuplot(ofstrOpt, "incrementalParam", incrementalParam, false);
    for (int i = 0; i < 4; i++) {
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
    drawGraphGnuplot("scripts/mggsa_test.gp", args);

#if defined(_MSC_VER)
    cin.get();
#endif

    return 0;
}
