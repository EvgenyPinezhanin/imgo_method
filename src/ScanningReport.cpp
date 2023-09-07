#include <reports/ScanningReport.h>

#include <algorithm>

#include "my_math.h" 

void ScanningReport::printOptProblem(std::ostream &stream, const OneDimensionalProblem &optProblem) const {
    stream << "[a; b] = [" << optProblem.getSearchArea().lowerBound << "; " <<
                              optProblem.getSearchArea().upBound << "]"<< "\n";
    vector<double> optimalPoints;
    optProblem.getOptimalPoints(optimalPoints);
    stream << "X* = (" << optimalPoints[0];
    int numberOptimalPoints = optimalPoints.size();
    for (int i = 1; i < numberOptimalPoints; i++) {
        stream << "; " << optimalPoints[i];
    }
    stream << ")\n";
    stream << "f(X*) = " << optProblem.getOptimalValue() << "\n";
}

void ScanningReport::printMethodParameters(
    std::ostream &stream,
    const opt::GeneralParametersNumericalOptMethod &parameters) const
{
    stream << "Method Parameters:" << "\n";
    stream << "Maximum of trials = " << parameters.maxTrials << "\n";
    stream << "Maximum of fevals = " << parameters.maxFevals << "\n";
    stream << "Accuracy = " << parameters.accuracy << "\n";
}

void ScanningReport::printResultMethod(std::ostream &stream, const ResultMethod &result) const {
    stream << "Result of method:" << "\n";
    stream << "Number of trials = " << result.numberTrials << "\n";
    stream << "Number of fevals = " << result.numberFevals << "\n";
    stream << "X = " << result.point << "\n";
    stream << "f(X) = " << result.value << "\n";
}

void ScanningReport::printErrorEstimate(
    std::ostream &stream,
    const OneDimensionalProblem &optProblem,
    const ResultMethod &result) const
{
    double point = result.point;
    vector<double> optimalPoints;
    optProblem.getOptimalPoints(optimalPoints);
    auto iter = std::min_element(optimalPoints.begin(), optimalPoints.end(),
        [&point] (const double &point1, const double &point2) {
            return std::abs(point1 - point) < std::abs(point2 - point);
        });
    stream << "|X* - X| = " << std::abs(*iter - point) << "\n";
    stream << "|f(X*) - f(X)| = " << std::abs(optProblem.getOptimalValue() - result.value) << "\n";
}
