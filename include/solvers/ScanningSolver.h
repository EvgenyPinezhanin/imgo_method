#ifndef SCANNING_SOLVER_H_
#define SCANNING_SOLVER_H_

#include <vector>

#include <solvers/IGeneralSolver.h>
#include <opt_methods/ScanningMethod.h>
#include <tasks/ScanningTask.h>
#include <result_methods/ResultMethod.h>
#include <reports/ScanningReport.h>
#include <trials/Trial.h>

class ScanningSolver:
    public IGeneralSolver<ScanningMethod, ScanningParameters, ScanningTask,
                          ResultMethod, ScanningReport, Trial, double> {
public:
    ScanningSolver():
        IGeneralSolver<ScanningMethod, ScanningParameters, ScanningTask,
                       ResultMethod, ScanningReport, Trial, double>() {};

};

#endif // SCANNING_SOLVER_H_
