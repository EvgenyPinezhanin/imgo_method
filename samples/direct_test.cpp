#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

#include <direct_method.h>
#include <tasks/task.h>
#include <output_results.h>

using namespace std;
/* 
const int functionNumber = 0; // 0 - f1, 1 - f2, ...
const int functionBlock = 0; // 0 - sample, 1 - test
const int displayType = 1; // 0 - application, 1 - png, 2 - png(notitle)

const int numberBlocks = 2;
const int numberFunctions[2] = { 4, 2 };
const vector<string> functionBlockName{ "sample", "test" };

double f1Sample(int n, const double *X, int *undefinedFlag, void *data) {
    vector<double> point{ X[0], X[1] };
    DataDirect *fData = static_cast<DataDirect*>(data);
    fData->numberFevals++;

    double f = 1.0 - X[0] - X[1];
    point.push_back(f);
    fData->points.push_back(point);

    return f;
}

double f2Sample(int n, const double *X, int *undefinedFlag, void *data) {
    vector<double> point{ X[0], X[1] };
    DataDirect *fData = static_cast<DataDirect*>(data);
    fData->numberFevals++;

    double f = (X[0] - 1.0) * (X[0] - 1.0) / 5.0 + (X[1] - 1.0) * (X[1] - 1.0) / 5.0;
    point.push_back(f);
    fData->points.push_back(point);

    return f;
}

double f3Sample(int n, const double *X, int *undefinedFlag, void *data) {
    vector<double> point{ X[0], X[1] };
    DataDirect *fData = static_cast<DataDirect*>(data);
    fData->numberFevals++;

    if (1.0 - X[0] - X[1] > 0.0) *undefinedFlag = 1;
    double f =  X[0] * X[0] / 5.0 + X[1] * X[1] / 5.0;
    point.push_back(f);
    fData->points.push_back(point);

    return f;
}

double f4Sample(int n, const double *X, int *undefinedFlag, void *data) {
    vector<double> point{ X[0], X[1] };
    DataDirect *fData = static_cast<DataDirect*>(data);
    fData->numberFevals++;

    if ((X[0] - 2.0) * (X[0] - 2.0) + (X[1] - 2.0) * (X[1] - 2.0) - 2.0 > 0.0) *undefinedFlag = 1;
    double f = X[0] * X[0] / 5.0 + X[1] * X[1] / 5.0;
    point.push_back(f);
    fData->points.push_back(point);

    return f;
}

double f1Test(int n, const double *X, int *undefinedFlag, void *data) {
    vector<double> point{ X[0], X[1] };
    DataDirect *fData = static_cast<DataDirect*>(data);
    fData->numberFevals++;

    bool constrained = true;
    vector<double> constraints(3);
    constraints[0] = 0.01 * (pow(X[0] - 2.2, 2) + pow(X[1] - 1.2, 2) - 2.25);
    constraints[1] = 100.0 * (1.0 - pow(X[0] - 2.0, 2) / 1.44 - pow(0.5 * X[1], 2));
    constraints[2] = 10.0 * (X[1] - 1.5 - 1.5 * sin(6.283 * (X[0] - 1.75)));
    for (int i = 0; i < 3; i++) {
        if (constraints[i] > 0.0) constrained = false;
    }
    if (!constrained) *undefinedFlag = 1;

    double f = -1.5 * X[0] * X[0] * exp(1.0 - X[0] * X[0] - 20.25 * pow((X[0] - X[1]), 2)) -
                pow(0.5 * (X[0] - 1.0) * (X[1] - 1.0), 4) * exp(2.0 - pow(0.5 * (X[0] - 1.0), 4) -
                pow(X[1] - 1.0, 4));
    point.push_back(f);
    fData->points.push_back(point);

    return f;
}

const double C[20] = { 75.1963666677, -3.8112755343, 0.1269366345, -0.0020567665, 0.000010345,
                      -6.8306567631, 0.0302344793, -0.0012813448, 0.0000352559, -0.0000002266,
                       0.2564581253, -0.0034604030, 0.0000135139, -28.1064434908, -0.0000052375,
                      -0.0000000063, 0.0000000007, 0.0003405462, -0.0000016638, -2.8673112392 };

double f2Test(int n, const double *X, int *undefinedFlag, void *data) {
    vector<double> point{ X[0], X[1] };
    DataDirect *fData = static_cast<DataDirect*>(data);
    fData->numberFevals++;

    bool constrained = true;
    vector<double> constraints(4);
    constraints[0] = 450.0 - X[0] * X[1];
    constraints[1] = (0.1 * X[0] - 1.0) * (0.1 * X[0] - 1.0) - X[1];
    constraints[2] = 8.0 * (X[0] - 40.0) - (X[1] - 30.0) * (X[1] - 55.0);
    constraints[3] = X[1] + (X[0] - 35.0) * (X[0] - 30.0) / 125.0 - 80.0;
    for (int i = 0; i < 4; i++) {
        if (constraints[i] > 0.0) constrained = false;
    }
    if (!constrained) *undefinedFlag = 1;

    double f = -(C[0] + C[1] * X[0] + C[2] * X[0] * X[0] + C[3] * pow(X[0], 3) + C[4] * pow(X[0], 4) + C[5] * X[1] +
                 C[6] * X[0] * X[1] + C[7] * X[0] * X[0] * X[1] + C[8] * pow(X[0], 3) * X[1] + C[9] * pow(X[0], 4) * X[1] +
                 C[10] * X[1] * X[1] + C[11] * pow(X[1], 3) + C[12] * pow(X[1], 4) + C[13] / (X[1] + 1) + C[14] * X[0] * X[0] * X[1] *X[1] +
                 C[15] * pow(X[0], 3) * X[1] * X[1] + C[16] * pow(X[0], 3) * pow(X[1], 3) + C[17] * X[0] * X[1] * X[1] +
                 C[18] * X[0] * pow(X[1], 3) + C[19] * exp(0.0005 * X[0] * X[1]));
    point.push_back(f);
    fData->points.push_back(point);

    return f;
} */

int main() {
/*    ofstream ofstr("output_data/direct_test.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstrOpt("output_data/direct_test_opt.txt");
    if (!ofstrOpt.is_open()) cerr << "File opening error\n";

    DataDirect fData;

    int n = 2;
    int maxFevals = 1000, maxIters = 1000;
    double magicEps = 1.0e-4;
    double methodAccuracy = 1.0e-5;
    double volumeReltol = methodAccuracy / sqrt(n);
    double sigmaReltol = 0.0;
    direct_algorithm algorithm = DIRECT_ORIGINAL;

    vector<TaskDirect> taskArray = {
        TaskDirect("f1(x, y) sample", f1Sample, &fData, n, vector<double>{ -4.0, -4.0 }, vector<double>{ 4.0, 4.0 },
                   vector<double>{ 4.0, 4.0 }, maxFevals, maxIters, magicEps, volumeReltol, sigmaReltol, nullptr, algorithm),
        TaskDirect("f2(x, y) sample", f2Sample, &fData, n, vector<double>{ -4.0, -4.0 }, vector<double>{ 4.0, 4.0 },
                   vector<double>{ 1.0, 1.0 }, maxFevals, maxIters, magicEps, volumeReltol, sigmaReltol, nullptr, algorithm),
        TaskDirect("f3(x, y) sample", f3Sample, &fData, n, vector<double>{ -1.0, -1.0 }, vector<double>{ 1.0, 1.0 },
                   vector<double>{ 0.5, 0.5 }, maxFevals, maxIters, magicEps, volumeReltol, sigmaReltol, nullptr, algorithm),
        TaskDirect("f4(x, y) sample", f4Sample, &fData, n, vector<double>{ 0.0, 0.0 }, vector<double>{ 3.0, 3.0 },
                   vector<double>{ 1.0, 1.0 }, maxFevals, maxIters, magicEps, volumeReltol, sigmaReltol, nullptr, algorithm),
        TaskDirect("f1(x, y) test", f1Test, &fData, n, vector<double>{ 0.0, -1.0 }, vector<double>{ 4.0, 3.0 },
                   vector<double>{ 0.942, 0.944 }, maxFevals, maxIters, magicEps, volumeReltol, sigmaReltol, nullptr, algorithm),
        TaskDirect("f2(x, y) test", f2Test, &fData, n, vector<double>{ 0.0, 0.0 }, vector<double>{ 80.0, 80.0 },
                   vector<double>{ 77.489, 63.858 }, maxFevals, maxIters, magicEps, volumeReltol, sigmaReltol, nullptr, algorithm)
    };

    DirectMethod direct;

    vector<double> X;
    double minF;
    int flagTmp;
    DataDirect *data, dataTmp;

    int taskArraySize = taskArray.size();
    for (int i = 0; i < taskArraySize; i++) {
        if (taskArray[i].used) {
            data = static_cast<DataDirect*>(taskArray[i].fData);
            data->numberFevals = 0;
            data->points.clear();

            direct.setF(taskArray[i].f);
            direct.setFData(taskArray[i].fData);
            direct.setN(taskArray[i].n);
            direct.setAB(taskArray[i].A, taskArray[i].B);
            direct.setMaxFevals(taskArray[i].maxFevals);
            direct.setMaxIters(taskArray[i].maxIters);
            direct.setMagicEps(taskArray[i].magicEps);
            direct.setVolumeReltol(taskArray[i].volumeReltol);
            direct.setSigmaReltol(taskArray[i].sigmaReltol);
            direct.setLogfile(taskArray[i].logfile);
            direct.setAlghorithm(taskArray[i].algorithm);

            direct.solve(X, minF);

            printResultDirect(taskArray[i].name, taskArray[i].n, taskArray[i].A, taskArray[i].B, taskArray[i].XOpt,
                              taskArray[i].f(taskArray[i].n, taskArray[i].XOpt.data(), &flagTmp, &dataTmp), taskArray[i].maxIters,
                              taskArray[i].maxFevals, taskArray[i].magicEps, taskArray[i].volumeReltol, taskArray[i].sigmaReltol,
                              taskArray[i].algorithm, data->numberFevals, X, minF);

            addPointGnuplot(ofstr, taskArray[i].XOpt, taskArray[i].f(taskArray[i].n, taskArray[i].XOpt.data(), &flagTmp, &dataTmp));
            addPointGnuplot(ofstr, X, taskArray[i].f(taskArray[i].n, X.data(), &flagTmp, &dataTmp));

            addPointsGnuplot(ofstr, data->points);
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
    drawGraphGnuplot("scripts/direct_test.gp", args);

#if defined( _MSC_VER )
    cin.get();
#endif
 */
    return 0;
}
