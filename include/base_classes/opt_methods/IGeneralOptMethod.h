#ifndef I_GENERAL_OPT_METHOD_H_
#define I_GENERAL_OPT_METHOD_H_

template <typename OptTaskType, typename ResultMethodType>
class IGeneralOptMethod {
protected:
    OptTaskType task;

    ResultMethodType result;

public:
    IGeneralOptMethod(const OptTaskType &_task, int _maxFevals)
                     : task(_task), result(ResultMethodType()) {};

    void setTask(const OptTaskType &_task) { task = _task; };
    OptTaskType getTask() const { return task; };

    virtual void solve(ResultMethodType &_result) = 0;
    virtual bool solveTest(ResultMethodType &_result) = 0;
};

#endif // I_GENERAL_OPT_METHOD_H_
