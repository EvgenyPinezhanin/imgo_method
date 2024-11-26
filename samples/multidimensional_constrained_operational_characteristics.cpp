#if defined( _MSC_VER )
    #define _CRT_SECURE_NO_WARNINGS
    #define PROC_BIND
#else
    #define PROC_BIND proc_bind(master)
#endif

#include <iostream>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <opt_methods/MggsaMethod.h>
#include <opt_methods/GsaMethod.h>
#include <test_opt_problems/FittingFamilyOptProblems.h>
#include <Solver.h>
#include <gnuplot/OutputFile.h>
#include <gnuplot/Script.h>
#include <MyMath.h>
#include <omp.h>

// #define CALC
#define DRAW

using OptMethod = GsaMethod<OneDimensionalSupportiveOptProblem>;
using MggsaParameters = MggsaMethod<FittingFamilyOptProblems<OptMethod>>::Parameters;
using GsaParameters = GsaMethod<OneDimensionalSupportiveOptProblem>::Parameters;
using TypeSolve = MggsaMethod<FittingFamilyOptProblems<OptMethod>>::TypeSolve;
using Task = opt::Task<FittingFamilyOptProblems<OptMethod>>;

const std::vector<std::string> methodNames{ "mggsa" };
const int familyNumber = 0; // 0 - Fitting
const int displayType = 0; // 0 - application, 1 - png, 2 - png(notitle)

void addOperationalCharacteristics(
    OutputFile &file,const std::vector<std::pair<size_t, double>> &operationalCharacteristics)
{
    size_t operationalCharacteristicsSize = operationalCharacteristics.size();
    for (size_t j = 0; j < operationalCharacteristicsSize; ++j) {
        file.addPoint(operationalCharacteristics[j].first, operationalCharacteristics[j].second, false);
    }
    file.close();
}

int main() {
    OutputFile varsFile("output_data/multidimensional_constrained_operational_characteristics/vars.txt");
    if (!varsFile.isOpen()) std::cerr << "vars.txt opening error\n";

    const size_t numberFamily = 1;

    double accuracyGsa = 0.01, reliabilityGsa = 2.5;
    size_t maxTrialsGsa = 10000, maxFevalsGsa = 10000;
    GsaParameters gsaParameters(accuracyGsa, 0.0, maxTrialsGsa, maxFevalsGsa, reliabilityGsa);

    GsaMethod<OneDimensionalSupportiveOptProblem> gsa;
    gsa.setParameters(gsaParameters);

    FittingFamilyOptProblems fittingFamilyOptProblems(gsa);

    size_t dimension = fittingFamilyOptProblems.getDimension();
    size_t familySize = fittingFamilyOptProblems.getFamilySize();

    std::vector<std::vector<std::vector<double>>> reliability{
        {
            std::vector<double>(dimension + 4, 10.0),
            std::vector<double>(dimension + 4, 7.0),
            std::vector<double>(dimension + 4, 5.0),
            std::vector<double>(dimension + 4, 3.0)
        }
    };
    // std::vector<size_t> key{ 3, 3, 3, 3 };
    std::vector<size_t> key{ 1, 1, 1, 1 };
    std::vector<double> d{ 0.01, 0.01, 0.01, 0.01 };

    std::vector<std::vector<size_t>> K{ { 0, 100000, 25 } };

    double error = 0.1;
    size_t maxFevals = 1000000;
    size_t density = 11, increment = 0;
    TypeSolve typeSolve = TypeSolve::SOLVE;
    MggsaParameters parameters(0.0, error, 0, maxFevals, std::vector<double>{},
                               0.0, density, 0, increment, typeSolve);

    MggsaMethod<FittingFamilyOptProblems<OptMethod>> mggsa;

    std::vector<Task> tasks{ Task( "FittingFamily", fittingFamilyOptProblems, parameters) };

    Solver<FittingFamilyOptProblems<OptMethod>> solver;
    std::vector<std::pair<size_t, double>> operationalCharacteristics;
    std::stringstream strReport, fileName;
    OutputFile operationalCharacteristicsFile;
    double workTime;

    // int chunk = 2;

    double totalStartTime = omp_get_wtime();
#if defined( CALC )
// #pragma omp parallel for schedule(dynamic, chunk) PROC_BIND num_threads(omp_get_num_procs()) collapse(2) \
//         shared(numberFamily, reliability, ) firstprivate(mggsa)
    for (size_t i = 0; i < numberFamily; ++i) {
        size_t reliabilitySize = reliability[i].size();

        for (size_t j = 0; j < reliabilitySize; ++j) {
            parameters.reliability = reliability[i][j];
            parameters.key = key[j];
            parameters.maxTrials = K[0][1];
            tasks[i].parameters = parameters;
            mggsa.setParameters(parameters);

            solver.calcOperationalCharacteristics(mggsa, tasks[i], K[0][0], K[0][1], K[0][2],
                                                  operationalCharacteristics, workTime, true);

            strReport << tasks[i].name << ", method: " << methodNames[0] << ", error = " << parameters.error
                      << ", max trials = " << parameters.maxTrials << ", max fevals = " << parameters.maxFevals
                      << ", reliability = " << parameters.reliability[0] << ", time: " << workTime << "\n";
            std::cout << strReport.str();
            strReport.str("");

            fileName << std::setprecision(2) << "output_data/multidimensional_constrained_operational_characteristics/"
                     << methodNames[0] << "/" << tasks[i].name << "_" << parameters.key << "_" << parameters.reliability[0];
            operationalCharacteristicsFile.open(fileName.str());
            if (!operationalCharacteristicsFile.isOpen()) std::cerr << fileName.str() << " opening error\n";
            fileName.str("");

            addOperationalCharacteristics(operationalCharacteristicsFile, operationalCharacteristics);
            operationalCharacteristicsFile.close();
        }
    }
#endif
    double totalEndTime = omp_get_wtime();
    std::cout << "Total time: " << totalEndTime - totalStartTime << std::endl;

    size_t sizeKey = key.size();
    varsFile.setVariable("numberKey", sizeKey, false);
    varsFile.initArray("methodNames", methodNames.size());
    varsFile.initArray("familyName", numberFamily);
    varsFile.initArray("r", sizeKey * numberFamily);
    varsFile.initArray("key", sizeKey);
    for (int i = 0; i < numberFamily; i++) {
        varsFile.setValueInArray("familyName", i + 1, tasks[i].name);
        for (int j = 0; j < reliability[i].size(); j++) {
            varsFile.setValueInArray("r", (i * sizeKey) + j + 1, reliability[i][j][0]);
        }
    }
    for (int i = 0; i < sizeKey; i++) {
        varsFile.setValueInArray("key", i + 1, key[i]);
    }
    for (int i = 0; i < methodNames.size(); i++) {
        varsFile.setValueInArray("methodNames", i + 1, methodNames[i]);
    }
    varsFile.close();

#if defined( DRAW )
    Script script("scripts/multidimensional_constrained_operational_characteristics.gp");
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
