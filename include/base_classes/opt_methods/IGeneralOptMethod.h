#ifndef I_GENERAL_OPT_METHOD_H_
#define I_GENERAL_OPT_METHOD_H_

namespace opt {
    template <typename OptProblemType, typename ParametersMethodType, typename ResultMethodType>
    class IGeneralOptMethod {
    protected:
        OptProblemType problem;
        ResultMethodType result;

    public:
        IGeneralOptMethod(const OptProblemType &_problem):
            problem(_problem),
            result(ResultMethodType())
        {};

        void setProblem(const OptProblemType &_problem) { problem = _problem; };
        OptProblemType getProblem() const { return problem; };

        virtual void setParameters(const ParametersMethodType &parameters) = 0;
        virtual void getParameters(ParametersMethodType &parameters) const = 0;

        virtual void solve(ResultMethodType &_result) = 0;
        virtual bool solveTest(ResultMethodType &_result) = 0;
    };
}

#endif // I_GENERAL_OPT_METHOD_H_
