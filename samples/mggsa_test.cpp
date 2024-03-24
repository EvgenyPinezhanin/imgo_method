#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <limits>

/* #include <opt_methods/MggsaMethod.h>
#include <gnuplot/OutputFile.h>
#include <gnuplot/Script.h>
#include <MyMath.h>

#define CALC
#define DRAW

const std::string methodName = "mggsa";
const int displayType = 1; // 0 - application, 1 - png, 2 - png(notitle)
const int problemBlock = 0; // 0 - sample, 1 - test
const int problemNumber = 2; // 0, 1, 2, ...*/

int main() {
/*     OutputFile trialsFile;

    double eps = 0.01, r = 2, d = 0.0;
    int n = 2, den = 10, key = 3;
    int maxTrials = 100000, maxFevals = 100000;

    std::vector<TaskMggsa> taskArray = { TaskMggsa(f1Sample, "f1(x, y) sample", n, 0, std::vector<double>{-4.0, -4.0},
                                                  std::vector<double>{4.0, 4.0}, std::vector<double>{4.0, 4.0},
                                                  std::vector<double>{}, eps, maxTrials, maxFevals, r, d, den, key, -1),
                                         TaskMggsa(f2Sample, "f2(x, y) sample", n, 0, std::vector<double>{-4.0, -4.0},
                                                   std::vector<double>{4.0, 4.0}, std::vector<double>{1.0, 1.0},
                                                   std::vector<double>{}, eps, maxTrials, maxFevals, r, d, den, key, -1),
                                         TaskMggsa(f3Sample, "f3(x, y) sample", n, 1, std::vector<double>{-1.0, -1.0},
                                                   std::vector<double>{1.0, 1.0}, std::vector<double>{0.5, 0.5},
                                                   std::vector<double>{}, eps, maxTrials, maxFevals, r, d, den, key, -1),
                                         TaskMggsa(f4Sample, "f4(x, y) sample", n, 1, std::vector<double>{0.0, 0.0},
                                                   std::vector<double>{3.0, 3.0}, std::vector<double>{1.0, 1.0},
                                                   std::vector<double>{}, eps, maxTrials, maxFevals, r, d, den, key, -1),

                                         TaskMggsa(f1Test, "f1(x, y) test", n, 3, std::vector<double>{ 0.0, -1.0 },
                                                   std::vector<double>{ 4.0, 3.0 }, std::vector<double>{ 0.942, 0.944 },
                                                   std::vector<double>{}, eps, maxTrials, maxFevals, r, d, den, key, -1),
                                         TaskMggsa(f2Test, "f2(x, y) test", n, 4, std::vector<double>{ 0.0, 0.0 },
                                                   std::vector<double>{ 80.0, 80.0 }, std::vector<double>{ 77.489, 63.858 },
                                                   std::vector<double>{}, eps, maxTrials, maxFevals, 3.3, 0.01, den, key, -1) };

    MggsaMethod mggsa;

    std::vector<double> X, L;
    std::vector<std::vector<double>> points;
    std::vector<opt::IndexTrial> trials;
    int numberTrials, numberFevals;

    int taskArraySize = taskArray.size();
    for (int i = 0; i < taskArraySize; i++) {
        if (taskArray[i].used) {
            if (i < 4) {
                trialsFile.open("output_data/" + methodName + "_test/" + functionBlockName[0] + "_" + std::to_string(i + 1));
            } else {
                trialsFile.open("output_data/" + methodName + "_test/" + functionBlockName[1] + "_" + std::to_string(i - 3));
            }
            if (!trialsFile.isOpen()) std::cerr << "trialsFile opening error\n";

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
                             taskArray[i].L, taskArray[i].X_opt, taskArray[i].f(taskArray[i].X_opt, taskArray[i].numberConstraints),
                             taskArray[i].maxTrials, taskArray[i].maxFevals, taskArray[i].eps, taskArray[i].r, taskArray[i].d,
                             taskArray[i].den, taskArray[i].key, taskArray[i].incr, numberTrials, numberFevals, L, X,
                             taskArray[i].f(X, taskArray[i].numberConstraints));

            trialsFile.addPoint(taskArray[i].X_opt, taskArray[i].f(taskArray[i].X_opt, taskArray[i].numberConstraints));
            trialsFile.addPoint(X, taskArray[i].f(X, taskArray[i].numberConstraints));

            mggsa.getPoints(points);
            trialsFile.addPoints(points);

            trialsFile.close();
        }
    }

#if defined( DRAW )
    Script script("scripts/multidimensional_constrained_test.gp");
    script.addArg(methodName);
    script.addArgs(std::vector<int>{ displayType, problemBlock, problemNumber });
    script.start();
    if (script.isError() == 2) std::cerr << "Error gnuplot\n";
    if (script.isError() == 1) std::cerr << "Error chmod\n";
#endif

#if defined( _MSC_VER )
    cin.get();
#endif */

    return 0;
}
