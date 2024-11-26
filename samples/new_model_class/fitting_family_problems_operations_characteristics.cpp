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

#include <random>

// #define CALC
#define DRAW

using OptMethod = GsaMethod<OneDimensionalSupportiveOptProblem>;
using MggsaParameters = MggsaMethod<FittingFamilyOptProblems<OptMethod>>::Parameters;
using GsaParameters = GsaMethod<OneDimensionalSupportiveOptProblem>::Parameters;
using TypeSolve = MggsaMethod<FittingFamilyOptProblems<OptMethod>>::TypeSolve;
using Task = opt::Task<FittingFamilyOptProblems<OptMethod>>;

const std::vector<std::string> methodNames{ "mggsa" };
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
    std::string rootDir = "output_data/new_model_class/fitting_family_problems_operational_characteristics";
    OutputFile varsFile(rootDir + "/vars.txt");
    if (!varsFile.isOpen()) std::cerr << "vars.txt opening error\n";

    double accuracyGsa = 0.01, reliabilityGsa = 2.5;
    size_t maxTrialsGsa = 10000, maxFevalsGsa = 10000;
    GsaParameters gsaParameters(accuracyGsa, 0.0, maxTrialsGsa, maxFevalsGsa, reliabilityGsa);

    GsaMethod<OneDimensionalSupportiveOptProblem> gsa;
    gsa.setParameters(gsaParameters);

    FittingFamilyOptProblems fittingFamilyOptProblems(gsa);

    size_t dimension = fittingFamilyOptProblems.getDimension();
    size_t numberConstraints = fittingFamilyOptProblems.getNumberConstraints();
    size_t familySize = fittingFamilyOptProblems.getFamilySize();

    size_t numberVariants = 3;
    std::vector<std::vector<double>> reliability {
        std::vector<double>(numberConstraints + 1, 3.5),
        std::vector<double>(numberConstraints + 1, 3.0),
        std::vector<double>(numberConstraints + 1, 2.0)
    };
    // std::vector<std::vector<double>> reliability {
    //     std::vector<double>(numberConstraints + 1, 4.0),
    //     std::vector<double>(numberConstraints + 1, 2.5),
    //     std::vector<double>(numberConstraints + 1, 5.0),
    //     std::vector<double>(numberConstraints + 1, 3.0)
    // };
    std::vector<size_t> key{ 3, 3, 3 };
    // std::vector<size_t> key{ 3, 3, 1, 1 };
    std::vector<double> d{ 0.01, 0.01, 0.01 };
    // std::vector<double> d{ 0.01, 0.01, 0.01, 0.01 };

    std::vector<std::vector<size_t>> K{ { 0, 20000, 25 } };

    double error = 0.1;
    size_t maxFevals = 1000000;
    size_t density = 11, increment = 0;
    TypeSolve typeSolve = TypeSolve::SOLVE;
    MggsaParameters parameters(0.0, error, 0, maxFevals, std::vector<double>{},
                               0.0, density, 0, increment, typeSolve);

    MggsaMethod<FittingFamilyOptProblems<OptMethod>> mggsa;

    Solver<FittingFamilyOptProblems<OptMethod>> solver;
    std::vector<std::pair<size_t, double>> operationalCharacteristics;
    std::stringstream strReport, fileName;
    OutputFile operationalCharacteristicsFile;
    double workTime;

    int chunk = 1;

    double totalStartTime = omp_get_wtime();
#if defined( CALC )
#pragma omp parallel for schedule(dynamic, chunk) PROC_BIND num_threads(omp_get_num_procs()) \
        shared(reliability, key, K) firstprivate(mggsa, numberVariants, parameters, fittingFamilyOptProblems)
    for (size_t i = 0; i < numberVariants; ++i) {
        Task task( "FittingFamily", fittingFamilyOptProblems, parameters);
        
        parameters.reliability = reliability[i];
        parameters.key = key[i];
        parameters.maxTrials = K[0][1];
        task.parameters = parameters;
        mggsa.setParameters(parameters);

        solver.calcOperationalCharacteristics(mggsa, task, K[0][0], K[0][1], K[0][2],
                                              operationalCharacteristics, workTime, true);

        strReport << task.name << ", method: " << methodNames[0] << ", error = " << parameters.error
                  << ", max trials = " << parameters.maxTrials << ", max fevals = " << parameters.maxFevals
                  << ", reliability = " << parameters.reliability[0] << ", time: " << workTime << "\n";
        std::cout << strReport.str();
        strReport.str("");

        fileName << std::setprecision(2) << rootDir << "/" << methodNames[0] << "/" << task.name << "_"
                 << parameters.key << "_" << parameters.reliability[0];
        operationalCharacteristicsFile.open(fileName.str());
        if (!operationalCharacteristicsFile.isOpen()) std::cerr << fileName.str() << " opening error\n";
        fileName.str("");

        addOperationalCharacteristics(operationalCharacteristicsFile, operationalCharacteristics);
        operationalCharacteristicsFile.close();
    }
#endif
    double totalEndTime = omp_get_wtime();
    std::cout << "Total time: " << totalEndTime - totalStartTime << std::endl;

    size_t numberVariantsDraw = 4;
    std::vector<std::vector<double>> reliabilityDraw {
        std::vector<double>(numberConstraints + 1, 4.0),
        // std::vector<double>(numberConstraints + 1, 3.5),
        // std::vector<double>(numberConstraints + 1, 3.0),
        // std::vector<double>(numberConstraints + 1, 2.5),
        std::vector<double>(numberConstraints + 1, 2.0),
        std::vector<double>(numberConstraints + 1, 5.0),
        std::vector<double>(numberConstraints + 1, 3.0)
    };
    std::vector<size_t> keyDraw{ 3, 3, /* 3, 3, 3, */ 1, 1 };

    varsFile.setVariable("numberKey", numberVariantsDraw, false);
    varsFile.setVariable("familyName", "FittingFamily");

    varsFile.initArray("methodNames", methodNames.size());
    varsFile.initArray("r", numberVariantsDraw);
    varsFile.initArray("key", numberVariantsDraw);
    for (int i = 0; i < numberVariantsDraw; ++i) {
        varsFile.setValueInArray("r", i + 1, reliabilityDraw[i][0]);
        varsFile.setValueInArray("key", i + 1, keyDraw[i]);
    }
    for (int i = 0; i < methodNames.size(); i++) {
        varsFile.setValueInArray("methodNames", i + 1, methodNames[i]);
    }
    varsFile.close();

#if defined( DRAW )
    Script script("scripts/new_model_class/fitting_family_problems_operational_characteristics.gp");
    script.addArgs(std::vector<int>{ displayType });
    script.start();
    if (script.isError() == 2) std::cerr << "Error gnuplot\n";
    if (script.isError() == 1) std::cerr << "Error chmod\n";
#endif

#if defined( _MSC_VER )
    cin.get();
#endif

	return 0;
}
