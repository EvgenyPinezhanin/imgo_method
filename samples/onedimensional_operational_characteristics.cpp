#include <iostream>
#include <string>
#include <vector>

#include <Hill/HillProblemFamily.hpp>
#include <Shekel/ShekelProblemFamily.hpp>
#include <opt_problems/OneDimensionalFamilyProblem.h>
#include <Solver.h>
#include <opt_methods/GsaMethod.h>
#include <opt_methods/PiyavskyMethod.h>
#include <opt_methods/ScanningMethod.h>
#include <gnuplot/OutputFile.h>
#include <gnuplot/Script.h>
#include <omp.h>

#define CALC
#define DRAW

using Task = opt::Task<OneDimensionalFamilyProblem>;
using Parameters = GsaMethod<OneDimensionalFamilyProblem>::Parameters;

const std::vector<std::string> methodNames{ "scanning", "piyavsky", "gsa" };
const int displayType = 1; // 0 - application, 1 - png, 2 - png(notitle)
const int familyNumber = 0; // 0 - Hill, 1 - Shekel

int main() {
    double error = 0.0001;
    int maxTrials = 100000, maxFevals = 100000;
    std::vector<std::vector<double>> reliability{ {3.0, 3.2, 3.4},
                                                  {3.1, 3.4, 3.7} };

    Parameters parameters(0.0, error, maxTrials, maxFevals, 0.0);

    const size_t numberFamily = 2;
    THillProblemFamily hillProblems;
    TShekelProblemFamily shekelProblems;

    std::vector<Task> problems{ Task(  "HillProblemFamily",   hillProblems, parameters),
                                Task("ShekelProblemFamily", shekelProblems, parameters) };

    size_t kStart = 0, kFinish = 500, kStep = 10;

    // Solver solver;

    double totalStartTime = omp_get_wtime();
#if defined( CALC )
    for (size_t i = 0; i < numberFamily; ++i) {
        
    }
#endif
    double totalEndTime = omp_get_wtime();
    std::cout << "Total time: " << totalEndTime - totalStartTime << std::endl;

#if defined( DRAW )
/*     VariablesFile variablesFile("output_data/onedimensional_operational_characteristics/vars.txt");
    if (!variablesFile.isOpen()) std::cerr << "Variables file opening error\n";

    
    initArrayGnuplot(ofstrOpt, "familyName", numberFamily);
    initArrayGnuplot(ofstrOpt, "r", r.size() * 3);
    for (int i = 0; i < numberFamily; i++) {
        setValueInArrayGnuplot(ofstrOpt, "familyName", i + 1, problems[i].shortName);
        for (int j = 0; j < r[i].size(); j++) {
            setValueInArrayGnuplot(ofstrOpt, "r", i * 3 + 1 + j, r[i][j]);
        }
    }

    variablesFile.closeFile(); */

    Script script("scripts/onedimensional_operational_characteristics.gp");
    script.addArgs(std::vector<int>{ displayType, familyNumber });
    script.start();
    if (script.isError() == 2) std::cerr << "Error gnuplot\n";
    if (script.isError() == 1) std::cerr << "Error chmod\n";
#endif

#if defined( _MSC_VER )
    cin.get();
#endif

	return 0;
}
