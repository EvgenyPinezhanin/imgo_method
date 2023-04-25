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
#include <task.h>
#include <output_results.h>

using namespace std;

const int functionNumber = 3; // 0 - f1, 1 - f2, 2 - f3, 3 - f4

double f1(double x) {
    return -4.0 * x + 1.0;
}

double f2(double x) {
    return 5.0 * x * x + 3.0 * x - 1.0;
}

double f3(double x) {
    return x * sin(x);
}

double f4(double x) {
    if (x != 0) {
        return x * sin(1 / x);
    } else {
        return 0.0;
    }
}

int main() {
    ofstream ofstr("output_data/gsa_sample.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";

    double eps = 0.0001, r = 2.0;
    int maxTrials = 100000, maxFevals = 100000;

    vector<TaskGsa> taskArray = { TaskGsa(f1, "f1(x) = -4.0 * x + 1.0", 3.0, 4.0, 4.0, 4.0, eps, maxTrials, maxFevals, r),
                                  TaskGsa(f2, "f2(x) = 5.0 * x * x + 3.0 * x - 1.0", -2.0, 2.0, -0.3, 23.0, eps, maxTrials, maxFevals, r),
                                  TaskGsa(f3, "f3(x) = x * sin(x)", 0.0, 20.0, 17.336, 18.955, eps , maxTrials, maxFevals, 2.1),
                                  TaskGsa(f4, "f4(x) = x * sin(1 / x)", -0.4, -0.05, -0.2225, 6.0 * M_PI, eps, maxTrials, maxFevals, r) };

    GsaMethod gsa;

    double x;
    int numberTrials, numberFevals;
    vector<Trial> trials;

    for (int i = 0; i < taskArray.size(); i++) {
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
                           taskArray[i].r, numberTrials, numberFevals, gsa.getLambda(), x, taskArray[i].f(x));

            addPointGnuplot(ofstr, taskArray[i].XOpt[0], taskArray[i].f(taskArray[i].XOpt[0]));
            addPointGnuplot(ofstr, x, taskArray[i].f(x));

            gsa.getTrialPoints(trials);
            addPointsGnuplot(ofstr, trials);
        }
    }
    ofstr.close();

    drawGraphGnuplot("scripts/gsa_sample.gp", functionNumber);

#if defined( _MSC_VER )
    cin.get();
#endif

    return 0;
}
