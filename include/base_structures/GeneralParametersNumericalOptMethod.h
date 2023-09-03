#ifndef GENERAL_PARAMETERS_NUMERICAL_OPT_METHOD_H_
#define GENERAL_PARAMETERS_NUMERICAL_OPT_METHOD_H_

namespace opt {
    struct GeneralParametersNumericalOptMethod {
        double accuracy, error;
        int maxTrials, maxFevals;

        GeneralParametersNumericalOptMethod(double _accuracy, double _error, int _maxTrials, int _maxFevals):
            accuracy(_accuracy),
            error(_error),
            maxTrials(_maxTrials),
            maxFevals(_maxFevals)
        {};
    };
}

#endif // GENERAL_PARAMETERS_NUMERICAL_OPT_METHOD_H_
