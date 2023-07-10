#ifndef MGGSA_RESULT_METHOD_H_
#define MGGSA_RESULT_METHOD_H_

#include <vector>

#include <base_structures/result_methods/GeneralNumericalResultMethod.h>

using std::vector;

enum class MggsaStopCriteria { accuracy, error, maxTrials, maxFevals, coincidePoints };

using MggsaResultMethod = GeneralNumericalResultMethod<MggsaStopCriteria, vector<double>>;

#endif // MGGSA_RESULT_METHOD_H_
