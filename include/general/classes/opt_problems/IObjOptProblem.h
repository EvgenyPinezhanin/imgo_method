#ifndef _I_OBJ_OPT_PROBLEM_H_
#define _I_OBJ_OPT_PROBLEM_H_

#include <vector>

namespace opt {
    template <typename ObjectType>
    class IObjOptProblem {
    protected:
        ObjectType object;

    public:
        IObjOptProblem(const ObjectType &_object) : object(_object) {};

        void getObject(ObjectType &_object) const { _object = object; };
    };
}

#endif // _I_OBJ_OPT_PROBLEM_H_
