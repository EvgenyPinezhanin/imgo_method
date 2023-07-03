#ifndef RESULT_METHOD_H_
#define RESULT_METHOD_H_

#include <vector>

using std::vector;

struct ResultMethod {
    int numberFevals;
    vector<double> X;

    ResultMethod(int _numberFevals = 0, vector<double> _X = vector<double>{})
                 : numberFevals(_numberFevals), X(_X) {};
};

#endif // RESULT_METHOD_H_
