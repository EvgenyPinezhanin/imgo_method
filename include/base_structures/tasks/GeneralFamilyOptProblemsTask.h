#ifndef GENERAL_FAMILY_OPT_PROBLEMS_TASK_H_
#define GENERAL_FAMILY_OPT_PROBLEMS_TASK_H_

#include <string>

using std::string;

template <typename FamilyOptProblemType>
struct GeneralFamilyOptProblemsTask {
    string name;

    FamilyOptProblemType familyOptProblem;

    bool use;

    GeneralFamilyOptProblemsTask(const string &_name, const OptFamilyOptProblemType &_familyOptProblem, bool _use)
                                : name(_name), familyOptProblem(_familyOptProblem), use(_use) {};
};

#endif // GENERAL_FAMILY_OPT_PROBLEMS_TASK_H_
