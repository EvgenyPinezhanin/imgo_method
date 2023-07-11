#ifndef PRINT_RESULT_H_
#define PRINT_RESULT_H_

#include <ostream>

#include <tasks/ScanningTask.h>
#include <result_methods/ResultMethod.h>

using std::ostream;

void printResult(ostream &ostr, const ScanningTask &task, const ResultMethod &result);

#endif // PRINT_RESULT_H_
