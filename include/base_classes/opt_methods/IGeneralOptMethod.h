#ifndef I_GENERAL_OPT_METHOD_H_
#define I_GENERAL_OPT_METHOD_H_

template <typename OptTaskType, typename ResultMethodType>
class IGeneralOptMethod {
protected:
    OptTaskType task;

    int maxFevals, numberFevals;

    virtual void setDataInResultMethod(ResultMethodType &result) const = 0;

public:
    IGeneralOptMethod(const OptTaskType &_task, int _maxFevals) : task(_task), maxFevals(_maxFevals), numberFevals(0) {};

    void setTask(const OptTaskType &_task) { task = _task; };
    OptTaskType getTask() const { return task; };

    void setMaxFevals(int _maxFevals) { maxFevals = _maxFevals; };
    int getMaxFevals() const { return maxFevals; };

    virtual void solve(ResultMethodType &result) = 0;
    virtual bool solveTest(ResultMethodType &result) = 0;
};

#endif // I_GENERAL_OPT_METHOD_H_
