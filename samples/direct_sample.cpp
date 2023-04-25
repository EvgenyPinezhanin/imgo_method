#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include <direct_method.h>
#include <task.h>
#include <output_results.h>

using namespace std;

const int functionNumber = 0; // 0 - f1, 1 - f2, 2 - f3, 3 - f4

double f1(int n, const double *X, int *undefinedFlag, void *data) {
    DataDirect *fData = static_cast<DataDirect*>(data);
    fData->numberFevals++;
    fData->points.push_back(vector<double>{X[0], X[1]});

    return 1.0 - X[0] - X[1];
}

double f2(int n, const double *X, int *undefinedFlag, void *data) {
    DataDirect *fData = static_cast<DataDirect*>(data);
    fData->numberFevals++;
    fData->points.push_back(vector<double>{X[0], X[1]});

    return (X[0] - 1.0) * (X[0] - 1.0) / 5.0 + (X[1] - 1.0) * (X[1] - 1.0) / 5.0;
}

double f3(int n, const double *X, int *undefinedFlag, void *data) {
    DataDirect *fData = static_cast<DataDirect*>(data);
    fData->numberFevals++;
    fData->points.push_back(vector<double>{X[0], X[1]});

    if (1.0 - X[0] - X[1] > 0.0) *undefinedFlag = 1;
    return X[0] * X[0] / 5.0 + X[1] * X[1] / 5.0;
}

double f4(int n, const double *X, int *undefinedFlag, void *data) {
    DataDirect *fData = static_cast<DataDirect*>(data);
    fData->numberFevals++;
    fData->points.push_back(vector<double>{X[0], X[1]});

    if ((X[0] - 2.0) * (X[0] - 2.0) + (X[1] - 2.0) * (X[1] - 2.0) - 2.0 > 0.0) *undefinedFlag = 1;
    return X[0] * X[0] / 5.0 + X[1] * X[1] / 5.0;
}

int main() {
    ofstream ofstr("output_data/direct_sample.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstrOpt("output_data/direct_sample_opt.txt");
    if (!ofstrOpt.is_open()) cerr << "File opening error\n";

    DataDirect fData;

    int n = 2;
    int maxFevals = 1000, maxIters = 1000;
    double magicEps = 1.0e-4;
    double methodAccuracy = 1.0e-5;
    double volumeReltol = methodAccuracy / sqrt(n);
    double sigmaReltol = 0.0;
    direct_algorithm algorithm = DIRECT_ORIGINAL;

    vector<TaskDirect> taskArray = { TaskDirect("f1", f1, &fData, n, vector<double>{ -4.0, -4.0 }, vector<double>{ 4.0, 4.0 },
                                                   vector<double>{ 4.0, 4.0 }, maxFevals, maxIters, magicEps, volumeReltol,
                                                   sigmaReltol, nullptr, algorithm),
                                     TaskDirect("f2", f2, &fData, n, vector<double>{ -4.0, -4.0 }, vector<double>{ 4.0, 4.0 },
                                                   vector<double>{ 1.0, 1.0 }, maxFevals, maxIters, magicEps, volumeReltol,
                                                   sigmaReltol, nullptr, algorithm),
                                     TaskDirect("f3", f3, &fData, n, vector<double>{ -1.0, -1.0 }, vector<double>{ 1.0, 1.0 },
                                                   vector<double>{ 0.5, 0.5 }, maxFevals, maxIters, magicEps, volumeReltol,
                                                   sigmaReltol, nullptr, algorithm),
                                     TaskDirect("f4", f4, &fData, n, vector<double>{ 0.0, 0.0 }, vector<double>{ 3.0, 3.0 },
                                                   vector<double>{ 1.0, 1.0 }, maxFevals, maxIters, magicEps, volumeReltol,
                                                   sigmaReltol, nullptr, algorithm) };

    DirectMethod direct;

    int flagTmp;
    DataDirect *data, dataTmp;
    vector<double> X;
    double minF;

    for (int i = 0; i < taskArray.size(); i++) {
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
                              taskArray[i].f(taskArray[i].n, taskArray[i].XOpt.data(), &flagTmp, &dataTmp), 
                              taskArray[i].maxIters, taskArray[i].maxFevals, taskArray[i].magicEps, taskArray[i].volumeReltol,
                              taskArray[i].sigmaReltol, taskArray[i].algorithm, data->numberFevals, X, minF);

            addPointGnuplot(ofstr, taskArray[i].XOpt, taskArray[i].f(taskArray[i].n, taskArray[i].XOpt.data(), &flagTmp, &dataTmp));
            addPointGnuplot(ofstr, X, taskArray[i].f(taskArray[i].n, X.data(), &flagTmp, &dataTmp));

            // addPointsGnuplot(ofstr, trials);
            for (int j = 0; j < data->points.size(); j++) {
                ofstr << data->points[j][0] << " " << data->points[j][1] << " " << 
                         taskArray[i].f(taskArray[i].n, data->points[j].data(), &flagTmp, &dataTmp) << endl;
            }
            ofstr << endl << endl;
        }
    }

    size_t sizeTaskArray = taskArray.size();
    string arraysName[] = { "AX","BX", "AY", "BY" };
    for (int i = 0; i < sizeTaskArray; i++) {
        initArrayGnuplot(ofstrOpt, arraysName[i], sizeTaskArray);
    }
    for (int i = 0; i < sizeTaskArray; i++) {
        for (int j = 0; j < 2; j++) {
            setValueInArrayGnuplot(ofstrOpt, arraysName[2 * j], i + 1, to_string(taskArray[i].A[j]));
            setValueInArrayGnuplot(ofstrOpt, arraysName[2 * j + 1], i + 1, to_string(taskArray[i].B[j]));
        }
    }
    ofstrOpt.close();

    drawGraphGnuplot("scripts/direct_sample.gp", functionNumber);

#if defined(_MSC_VER)
    cin.get();
#endif

    return 0;
}
