#include <iostream>
#include <utility>
#include <string>
#include <sstream>
#include <vector>

#include <Hill/HillProblemFamily.hpp>
#include <Shekel/ShekelProblemFamily.hpp>
#include <opt_problems/OneDimensionalFamilyOptProblemsGCGen.h>
#include <opt_methods/GsaMethod.h>
#include <opt_methods/PiyavskyMethod.h>
#include <opt_methods/ScanningMethod.h>
#include <Solver.h>
#include <gnuplot/OutputFile.h>
#include <gnuplot/Script.h>
#include <omp.h>

#define CALC
#define DRAW

using FamilyOptProblems = OneDimensionalFamilyOptProblemsGCGen;
using Task = opt::Task<FamilyOptProblems>;
using Parameters = GsaMethod<FamilyOptProblems>::Parameters;
using OptMethod = opt::IGeneralNumericalOptMethod<GsaMethod<FamilyOptProblems>::Trial, FamilyOptProblems>;

const std::vector<std::string> methodNames{ "scanning", "piyavsky", "gsa" };
const int displayType = 0; // 0 - application, 1 - png, 2 - png(notitle)
const int familyNumber = 0; // 0 - Hill, 1 - Shekel

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
    size_t kStart = 0, kFinish = 500, kStep = 10;

    double error = 0.0001;
    int maxTrials = kFinish, maxFevals = kFinish;
    std::vector<std::vector<double>> reliability{ {3.0, 3.2, 3.4},
                                                  {3.1, 3.4, 3.7} };

    Parameters parameters(0.0, error, maxTrials, maxFevals, 0.0);

    const size_t numberMethods = 3;
    ScanningMethod<FamilyOptProblems> scanning;
    PiyavskyMethod<FamilyOptProblems> piyavsky;
    GsaMethod<FamilyOptProblems> gsa;

    std::vector<OptMethod*> methods{ &scanning, &piyavsky, &gsa };

    const size_t numberFamily = 2;
    THillProblemFamily hillProblems;
    TShekelProblemFamily shekelProblems;

    std::vector<Task> problems{ Task(  "HillFamily",   &hillProblems, parameters),
                                Task("ShekelFamily", &shekelProblems, parameters) };

    Solver<FamilyOptProblems> solver;
    std::vector<std::pair<size_t, double>> operationalCharacteristics;
    size_t operationalCharacteristicsSize, reliabilitySize;
    std::stringstream strReport, fileName;
    OutputFile operationalCharacteristicsFile;
    double workTime;

    double totalStartTime = omp_get_wtime();
#if defined( CALC )
    for (size_t i = 0; i < numberFamily; ++i) {
        solver.calcOperationalCharacteristics(*methods[0], problems[i], kStart, kFinish, kStep,
                                              operationalCharacteristics, workTime);

        strReport << problems[i].name << ", method: " << methodNames[0] << ", error = " << parameters.error
                  << ", max trials = " << parameters.maxTrials << ", max fevals = " << parameters.maxFevals
                  << ", time: " << workTime << "\n";
        std::cout << strReport.str();
        strReport.str("");

        fileName << std::setprecision(2) << "output_data/onedimensional_operational_characteristics/" + methodNames[0] + "/" + problems[i].name;
        operationalCharacteristicsFile.open(fileName.str());
        if (!operationalCharacteristicsFile.isOpen()) std::cerr << fileName.str() << " opening error\n";
        fileName.str("");

        addOperationalCharacteristics(operationalCharacteristicsFile, operationalCharacteristics);
        operationalCharacteristicsFile.close();

        for (size_t j = 1; j < numberMethods; ++j) {
            reliabilitySize = reliability[j - 1].size();
            for (size_t k = 0; k < reliabilitySize; ++k) {
                parameters.reliability = reliability[j - 1][k];
                methods[j]->setParameters(parameters);

                solver.calcOperationalCharacteristics(*methods[j], problems[i], kStart, kFinish, kStep,
                                                      operationalCharacteristics, workTime);

                strReport << problems[i].name << ", method: " << methodNames[j] << ", error = " << parameters.error
                          << ", max trials = " << parameters.maxTrials << ", max fevals = " << parameters.maxFevals
                          << ", reliability = " << parameters.reliability << ", time: " << workTime << "\n";
                std::cout << strReport.str();
                strReport.str("");

                fileName << std::setprecision(2) << "output_data/onedimensional_operational_characteristics/"
                         << methodNames[j] << "/" << problems[i].name << "_" << parameters.reliability;
                operationalCharacteristicsFile.open(fileName.str());
                if (!operationalCharacteristicsFile.isOpen()) std::cerr << fileName.str() << " opening error\n";
                fileName.str("");

                addOperationalCharacteristics(operationalCharacteristicsFile, operationalCharacteristics);
                operationalCharacteristicsFile.close();
            }
        }
    }
#endif
    double totalEndTime = omp_get_wtime();
    std::cout << "Total time: " << totalEndTime - totalStartTime << std::endl;

#if defined( DRAW )
    fileName << "output_data/onedimensional_operational_characteristics/vars.txt";
    OutputFile variablesFile(fileName.str());
    if (!variablesFile.isOpen()) std::cerr << fileName.str() << " file opening error\n";

    variablesFile.initArray("methodNames", methodNames);
    variablesFile.initArray("familyName", numberFamily);
    variablesFile.initArray("r", reliability.size() * 3);
    for (int i = 0; i < numberFamily; i++) {
        variablesFile.setValueInArray("familyName", i + 1, problems[i].name);
        for (int j = 0; j < reliability[i].size(); j++) {
            variablesFile.setValueInArray("r", i * 3 + 1 + j, reliability[i][j]);
        }
    }

    variablesFile.close();

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
