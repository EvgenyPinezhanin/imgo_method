#ifndef STRONGIN_RESULT_METHOD_H_
#define STRONGIN_RESULT_METHOD_H_

#include <ResultMethod.h>

enum class StopCriteria { accuracy, maxTrials, maxFevals, coincidePoints };

struct StronginResultMethod : public ResultMethod {
    int numberTrials;
    StopCriteria stop;

    StronginResultMethod(int _numberTrials = 0, int _numberFevals = 0, StopCriteria _stop = StopCriteria::accuracy,
                         vector<double> _X = vector<double>{})
                         : ResultMethod(_numberFevals, _X), numberTrials(_numberTrials), stop(_stop) {};
};

#endif // STRONGIN_RESULT_METHOD_H_
