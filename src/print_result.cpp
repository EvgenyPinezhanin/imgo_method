#include <print_result.h>

#include <iomanip>

#include <my_math.h>

using std::setprecision;
using std::endl;
using std::abs;

void printResult(ostream &ostr, const ScanningTask &task, const ResultMethod &result) {
    const auto defaultPrecision = ostr.precision();
    ostr << setprecision(8);

    ostr << "Task name: " << task.name << "\n";
    ostr << "[a; b] = [" << task.optTask.getSearchArea().getLowerBound() << "; " <<
                            task.optTask.getSearchArea().getUpBound() << "]"<< "\n";
    double xOpt = task.optTask.getOptimalPoint();
    ostr << "X* = " << xOpt << "\n";
    double fXOpt = task.optTask.computeObjFunction(xOpt);
    ostr << "f(X*) = " << fXOpt << "\n";

    ostr << "Method Parameters:" << "\n";
    ostr << "Maximum of trials = " << task.maxTrials << "\n";
    ostr << "Maximum of fevals = " << task.maxFevals << "\n";
    ostr << "Accuracy = " << task.accuracy << "\n";

    ostr << "Result of method:" << "\n";
    ostr << "Number of trials = " << result.numberTrials << "\n";
    ostr << "Number of fevals = " << result.numberFevals << "\n";
    ostr << "X = " << result.solution << "\n";
    double fX = task.optTask.computeObjFunction(result.solution);
    ostr << "f(X) = " << fX << "\n";
    ostr << "|X* - X| = " << abs(xOpt - result.solution) << "\n";
    ostr << "|f(X*) - f(X)| = " << abs(fXOpt - fX) << "\n";
    ostr << endl;

    ostr << setprecision(defaultPrecision);
}
