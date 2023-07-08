#ifndef GENERAL_RESULT_METHOD_H_
#define GENERAL_RESULT_METHOD_H_

template <typename SolutionType>
struct GeneralResultMethod {
    int numberFevals;
    SolutionType solution;

    GeneralResultMethod(int _numberFevals, SolutionType _solution)
                        : numberFevals(_numberFevals), solution(_solution) {};
};

#endif // GENERAL_RESULT_METHOD_H_
