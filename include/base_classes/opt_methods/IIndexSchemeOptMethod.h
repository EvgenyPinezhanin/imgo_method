#ifndef I_INDEX_SCHEME_OPT_METHOD_H_
#define I_INDEX_SCHEME_OPT_METHOD_H_

#include <vector>

#include <base_structures/trials/IndexTrial.h>

using std::vector;

namespace opt {
    class IIndexSchemeOptMethod {
    protected:
        vector<double> reliability;
        vector<double> constantsEstimation;

        vector<vector<IndexTrial>> I;
        vector<double> zValues;

        virtual void estimatingConstants() = 0;

        virtual void calcI() = 0;
        virtual void calcZValues() = 0;

    public:
        IIndexSchemeOptMethod(const vector<double> &_reliability)
                              : reliability(_reliability), constantsEstimation(0), I(0), zValues(0) {};

        void setReliability(const vector<double> &_reliability) { reliability = _reliability; };
        void getReliability(vector<double>& _reliability) const { _reliability = reliability; };

        void getConstantEstimation(vector<double> &_constantsEstimation) const {
            _constantsEstimation = constantsEstimation;
        };
    };
}

#endif // I_INDEX_SCHEME_OPT_METHOD_H_
