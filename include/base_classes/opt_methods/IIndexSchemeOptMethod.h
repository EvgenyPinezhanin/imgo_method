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
        virtual void calcZValues() = 0;

    public:
        IIndexSchemeOptMethod(const vector<double> &_reliability)
            : reliability(_reliability), constantsEstimation(reliability.size()),
              I(reliability.size()), zValues(reliability.size()) {};

        void setReliability(const vector<double> &_reliability) { 
            reliability = _reliability;
            constantsEstimation.resize(reliability.size());
            I.resize(reliability.size());
            zValues.resize(reliability.size());
        };
        void getReliability(vector<double>& _reliability) const { _reliability = reliability; };

        void getConstantEstimation(vector<double> &_constantsEstimation) const {
            _constantsEstimation = constantsEstimation;
        };
    };
}

#endif // I_INDEX_SCHEME_OPT_METHOD_H_
