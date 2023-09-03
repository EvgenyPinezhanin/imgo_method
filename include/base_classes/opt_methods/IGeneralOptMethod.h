#ifndef I_GENERAL_OPT_METHOD_H_
#define I_GENERAL_OPT_METHOD_H_

namespace opt {
    template <typename OptProblemType, typename ParametersMethod, typename ResultMethodType>
    class IGeneralOptMethod {
    protected:
        OptProblemType problem;
        ParametersMethod parameters;
        ResultMethodType result;

    public:
        IGeneralOptMethod(const OptProblemType &_problem, const ParametersMethod &_parameters):
            problem(_problem),
            parameters(_parameters),
            result(ResultMethodType())
        {};

        void setProblem(const OptProblemType &_problem) { problem = _problem; };
        OptProblemType getProblem() const { return problem; };

        void setParameters(const ParametersMethod &_parameters) { parameters = _parameters; };
        ParametersMethod getParameters() const { return parameters; };

        virtual void solve(ResultMethodType &_result) = 0;
        virtual bool solveTest(ResultMethodType &_result) = 0;
    };
}

#endif // I_GENERAL_OPT_METHOD_H_
