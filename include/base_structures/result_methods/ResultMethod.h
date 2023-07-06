#ifndef RESULT_METHOD_H_
#define RESULT_METHOD_H_

template <typename SolutionType>
struct ResultMethod {
    int numberFevals;
    SolutionType solution;

    ResultMethod(int _numberFevals, SolutionType _solution) : numberFevals(_numberFevals), solution(_solution) {};
};

#endif // RESULT_METHOD_H_
