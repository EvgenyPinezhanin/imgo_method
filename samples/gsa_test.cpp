#include <iostream>
#include <fstream>
#include <vector>

#if defined( _MSC_VER )
    #define _USE_MATH_DEFINES
    #include <math.h>
#else
    #include <cmath>
#endif

#include <gsa.h>
// #include <task.h>
// #include <output_results.h>

using namespace std;

const int functionNumber = 0; // 0 - f1, 1 - f2, ...
const int functionBlock = 0; // 0 - sample, 1 - test
const int displayType = 1; // 0 - application, 1 - png, 2 - png(notitle)

const int numberBlocks = 2;
const int numberFunctions[2] = { 4, 20 };
const vector<string> functionBlockName{ "sample", "test" };

double f1Sample(double x) {
    return -4.0 * x + 1.0;
}

double f2Sample(double x) {
    return 5.0 * x * x + 3.0 * x - 1.0;
}

double f3Sample(double x) {
    return x * sin(x);
}

double f4Sample(double x) {
    return  x != 0 ? x * sin(1 / x) : 0.0;
}

double f1Test(double x) {
    return -(-1.0 / 6.0 * pow(x, 6) + 52.0 / 25.0 * pow(x, 5) - 39.0 / 80.0 * pow(x, 4) - 
             71.0 / 10.0 * pow(x, 3) + 79.0 / 20.0 * x * x + x - 1.0 / 10.0);
}

double f2Test(double x) {
    return -(-sin(x) - sin(10.0 / 3.0 * x));
}

double f3Test(double x) {
    double sum = 0;
    for (int i = 1; i <= 5; i++) {
        sum += i * sin((i + 1.0) * x  + i);
    }
    return -sum;
}

double f4Test(double x) {
    return -(16.0 * x * x - 24.0 * x + 5.0) * exp(-x);
}

double f5Test(double x) {
    return -(-3.0 * x + 1.4) * sin(18.0 * x);
}

double f6Test(double x) {
    return -((x + sin(x)) * exp(-x * x));
}

double f7Test(double x) {
    return -(-sin(x) - sin(10.0 / 3.0 * x) - log(x) + 0.84 * x - 3.0);
}

double f8Test(double x) {
    double sum = 0;
    for (int i = 1; i <= 5; i++) {
        sum += i * cos((i + 1.0) * x  + i);
    }
    return -sum;
}

double f9Test(double x) {
    return -(-sin(x) - sin(2.0 / 3.0 * x));
}

double f10Test(double x) {
    return -(x * sin(x));
}

double f11Test(double x) {
    return -(-2.0 * cos(x) - cos(2.0 * x));
}

double f12Test(double x) {
    return -(-pow(sin(x), 3) - pow(cos(x), 3));
}

double f13Test(double x) {
    return x * x - 1 < 0 ? -(pow(x, 2.0 / 3.0) + pow(-(x * x - 1), 1.0 / 3.0)) :
                           -(pow(x, 2.0 / 3.0) - pow(x * x - 1, 1.0 / 3.0));
}

double f14Test(double x) {
    return -(exp(-x) * sin(2 * M_PI * x));
}

double f15Test(double x) {
    return -((-x * x + 5.0 * x - 6.0) / (x * x + 1));
}

double f16Test(double x) {
    return -(-2.0 * (x - 3) * (x - 3) - exp(- x * x / 2));
}

double f17Test(double x) {
    return -(-pow(x, 6) + 15.0 * pow(x, 4) - 27.0 * x * x - 250.0);
}

double f18Test(double x) {
    return x <= 3.0 ? (x - 2.0) * (x - 2.0) : -(-2.0 * log(x - 2.0) - 1.0);
}

double f19Test(double x) {
    return -(x - sin(3.0 * x) + 1.0);
}

double f20Test(double x) {
    return -(x - sin(x)) * exp(- x * x);
}

int main() {
/*     ofstream ofstr("output_data/gsa_test.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstrOpt("output_data/gsa_test_opt.txt");
    if (!ofstrOpt.is_open()) cerr << "File opening error\n";

    double eps = 0.001, r = 2.0;
    int maxTrials = 100000, maxFevals = 100000;

    vector<TaskGsa> taskArray = { TaskGsa(f1Sample, "f1(x) sample", 3.0, 4.0, 4.0, 4.0, eps, maxTrials, maxFevals, r),
                                  TaskGsa(f2Sample, "f2(x) sample", -2.0, 2.0, -0.3, 23.0, eps, maxTrials, maxFevals, r),
                                  TaskGsa(f3Sample, "f3(x) sample", 0.0, 20.0, 17.336, 18.955, eps , maxTrials, maxFevals, 2.1),
                                  TaskGsa(f4Sample, "f4(x) sample", -0.4, -0.05, -0.2225, 6.0 * M_PI, eps, maxTrials, maxFevals, r),
                                  TaskGsa(f1Test, "f1(x) test", -1.5, 11.0, 10.0, 13870.0, eps, maxTrials, maxFevals, r),
                                  TaskGsa(f2Test, "f2(x) test", 2.7, 7.5, 5.145735, 4.29, eps, maxTrials, maxFevals, r),
                                  TaskGsa(f3Test, "f3(x) test", -10.0, 10.0, -0.49139, 67.0, eps, maxTrials, maxFevals, r),
                                  TaskGsa(f4Test, "f4(x) test", 1.9, 3.9, 2.868, 3.0, eps, maxTrials, maxFevals, r),
                                  TaskGsa(f5Test, "f5(x) test", 0.0, 1.2, 0.96609, 36.0, eps, maxTrials, maxFevals, r),
                                  TaskGsa(f6Test, "f6(x) test", -10.0, 10.0, 0.67956, 2.5, eps, maxTrials, maxFevals, r),
                                  TaskGsa(f7Test, "f7(x) test", 2.7, 7.5, 5.19978, 6.0, eps, maxTrials, maxFevals, r),
                                  TaskGsa(f8Test, "f8(x) test", -10.0, 10.0, 5.48286, 67.0, eps, maxTrials, maxFevals, r),
                                  TaskGsa(f9Test, "f9(x) test", 3.1, 20.4, 17.039, 1.7, eps, maxTrials, maxFevals, r),
                                  TaskGsa(f10Test, "f10(x) test", 0.0, 10.0, 7.9787, 11.0, eps, maxTrials, maxFevals, r),
                                  TaskGsa(f11Test, "f11(x) test", -1.57, 6.28, 2.094, 3.0, eps, maxTrials, maxFevals, r),
                                  TaskGsa(f12Test, "f12(x) test", 0.0, 6.28, 4.712, 2.2, eps, maxTrials, maxFevals, r),
                                  TaskGsa(f13Test, "f13(x) test", 0.001, 0.99, 0.7071, 8.5, eps, maxTrials, maxFevals, r),
                                  TaskGsa(f14Test, "f14(x) test", 0.0, 4.0, 0.224885, 6.5, eps, maxTrials, maxFevals, r),
                                  TaskGsa(f15Test, "f15(x) test", -5.0, 5.0, 2.4142, 6.5, eps, maxTrials, maxFevals, r),
                                  TaskGsa(f16Test, "f16(x) test", -3.0, 3.0, 3.0, 85.0, eps, maxTrials, maxFevals, r),
                                  TaskGsa(f17Test, "f17(x) test", -4.0, 4.0, -3.0, 2520.0, eps, maxTrials, maxFevals, r),
                                  TaskGsa(f18Test, "f18(x) test", 0.0, 6.0, 2.0, 4.0, eps, maxTrials, maxFevals, r),
                                  TaskGsa(f19Test, "f19(x) test", 0.0, 6.5, 5.87287, 4.0, eps, maxTrials, maxFevals, r),
                                  TaskGsa(f20Test, "f20(x) test", -10.0, 10.0, 1.195137, 1.3, eps, maxTrials, maxFevals, r) };

    GsaMethod gsa;

    double x;
    int numberTrials, numberFevals;
    vector<Trial> trials;

    int taskArraySize = taskArray.size();
    for (int i = 0; i < taskArraySize; i++) {
        if (taskArray[i].used) {
            gsa.setF(taskArray[i].f);
            gsa.setAB(taskArray[i].A[0], taskArray[i].B[0]);
            gsa.setEps(taskArray[i].eps);
            gsa.setMaxTrials(taskArray[i].maxTrials);
            gsa.setMaxFevals(taskArray[i].maxFevals);
            gsa.setR(taskArray[i].r);

            gsa.solve(numberTrials, numberFevals, x);

            printResultGsa(taskArray[i].name, taskArray[i].A[0], taskArray[i].B[0], taskArray[i].L[0], taskArray[i].XOpt[0],
                           taskArray[i].f(taskArray[i].XOpt[0]), taskArray[i].maxTrials, taskArray[i].maxFevals, taskArray[i].eps,
                           taskArray[i].r, numberTrials, numberFevals, gsa.getL(), x, taskArray[i].f(x));

            addPointGnuplot(ofstr, taskArray[i].XOpt[0], taskArray[i].f(taskArray[i].XOpt[0]));
            addPointGnuplot(ofstr, x, taskArray[i].f(x));

            gsa.getTrialPoints(trials);
            addPointsGnuplot(ofstr, trials);
        }
    }
    ofstr.close();

    initArrayGnuplot(ofstrOpt, "A", taskArraySize);
    initArrayGnuplot(ofstrOpt, "B", taskArraySize);
    for (int i = 0; i < taskArraySize; i++) {
        setValueInArrayGnuplot(ofstrOpt, "A", i + 1, taskArray[i].A[0], false);
        setValueInArrayGnuplot(ofstrOpt, "B", i + 1, taskArray[i].B[0], false);
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
    drawGraphGnuplot("scripts/gsa_test.gp", args);

#if defined( _MSC_VER )
    cin.get();
#endif
 */
    return 0;
}
