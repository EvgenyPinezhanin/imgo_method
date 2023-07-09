#ifndef GENERAL_RESULT_METHOD_H_
#define GENERAL_RESULT_METHOD_H_

template <typename SolutionType>
struct GeneralResultMethod {
    SolutionType solution;

    GeneralResultMethod(SolutionType _solution = SolutionType()) : solution(_solution) {};
};

#endif // GENERAL_RESULT_METHOD_H_
