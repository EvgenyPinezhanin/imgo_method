#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include <mggsa.h>
#include <task.h>
#include <output_results.h>

using namespace std;

const int functionNumber = 0; // 0 - f1, 1 - f2, ...
const int functionBlock = 0; // 0 - sample, 1 - test
const int displayType = 1; // 0 - application, 1 - png, 2 - png(notitle)

const int numberBlocks = 2;
const int numberFunctions[2] = { 4, 2 };
const vector<string> functionBlockName{ "sample", "test" };

double f1Sample(vector<double> x, int j) {
    switch (j) {
        case 1: return 1.0 - x[0] - x[1];
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f2Sample(vector<double> x, int j) {
    switch (j) {
        case 1: return (x[0] - 1.0) * (x[0] - 1.0) / 5.0 + (x[1] - 1.0) * (x[1] - 1.0) / 5.0;
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f3Sample(vector<double> x, int j) {
    switch (j) {
        case 1: return 1.0 - x[0] - x[1];
        case 2: return x[0] * x[0] / 5.0 + x[1] * x[1] / 5.0;
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f4Sample(vector<double> x, int j) {
    switch (j) {
        case 1: return (x[0] - 2.0) * (x[0] - 2.0) + (x[1] - 2.0) * (x[1] - 2.0) - 2.0;
        case 2: return x[0] * x[0] / 5.0 + x[1] * x[1] / 5.0;
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f1Test(vector<double> x, int j) {
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

double f2Test(vector<double> x, int j) {
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

    double eps = 0.01, r = 1.7, d = 0.005;
    int n = 2, den = 8, key = 3;
    int maxTrials = 100000, maxFevals = 100000;

    vector<TaskMggsa> taskArray = { TaskMggsa(f1Sample, "f1(x, y) sample", n, 0, vector<double>{-4.0, -4.0}, vector<double>{4.0, 4.0},
                                              vector<double>{4.0, 4.0}, vector<double>{}, eps, maxTrials, maxFevals, r, d, den, key, -1),
                                    TaskMggsa(f2Sample, "f2(x, y) sample", n, 0, vector<double>{-4.0, -4.0}, vector<double>{4.0, 4.0},
                                              vector<double>{1.0, 1.0}, vector<double>{}, eps, maxTrials, maxFevals, r, d, den, key, -1),
                                    TaskMggsa(f3Sample, "f3(x, y) sample", n, 1, vector<double>{-1.0, -1.0}, vector<double>{1.0, 1.0},
                                              vector<double>{0.5, 0.5}, vector<double>{}, eps, maxTrials, maxFevals, r, d, den, key, -1),
                                    TaskMggsa(f4Sample, "f4(x, y) sample", n, 1, vector<double>{0.0, 0.0}, vector<double>{3.0, 3.0},
                                              vector<double>{1.0, 1.0}, vector<double>{}, eps, maxTrials, maxFevals, r, d, den, key, -1),
                                    TaskMggsa(f1Test, "f1(x, y) test", n, 3, vector<double>{ 0.0, -1.0 }, vector<double>{ 4.0, 3.0 },
                                              vector<double>{ 0.942, 0.944 }, vector<double>{}, eps, maxTrials, maxFevals, r, d, den,
                                              key, -1),
                                    TaskMggsa(f2Test, "f2(x, y) test", n, 4, vector<double>{ 0.0, 0.0 }, vector<double>{ 80.0, 80.0 },
                                              vector<double>{ 77.489, 63.858 }, vector<double>{}, eps, maxTrials, maxFevals, 3.3, 0.01,
                                              den, key, -1) };

    MggsaMethod mggsa;

    vector<double> X, L;
    vector<vector<double>> points;
    vector<TrialConstrained> trials;
    int numberTrials, numberFevals;

    int taskArraySize = taskArray.size();
    for (int i = 0; i < taskArraySize; i++) {
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
            addPointsGnuplot(ofstr, points);
        }
    }
    ofstr.close();

    initArrayGnuplot(ofstrOpt, "A", 2 * taskArraySize);
    initArrayGnuplot(ofstrOpt, "B", 2 * taskArraySize);
    for (int i = 0; i < taskArraySize; i++) {
        for (int j = 0; j < 2; j++) {
            setValueInArrayGnuplot(ofstrOpt, "A", 2 * i + 1 + j, taskArray[i].A[j], false);
            setValueInArrayGnuplot(ofstrOpt, "B", 2 * i + 1 + j, taskArray[i].B[j], false);
        }
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
    drawGraphGnuplot("scripts/mggsa_test.gp", args);

#if defined(_MSC_VER)
    cin.get();
#endif

    return 0;
}
