#ifndef _I_INDEX_SCHEME_OPT_METHOD_H_
#define _I_INDEX_SCHEME_OPT_METHOD_H_

#include <vector>

#include <general/structures/trials/IndexTrial.h>

namespace opt {
    class IIndexSchemeOptMethod {
    protected:
        std::vector<double> reliability;
        std::vector<double> constantsEstimation;

        std::vector<std::vector<IndexTrial>> I;
        std::vector<double> zValues;

        virtual void estimatingConstants() = 0;
        virtual void calcZValues() = 0;

    public:
        IIndexSchemeOptMethod(const std::vector<double> &_reliability)
            : reliability(_reliability), constantsEstimation(reliability.size()),
              I(reliability.size()), zValues(reliability.size()) {};

        void setReliability(const std::vector<double> &_reliability) { 
            reliability = _reliability;
            constantsEstimation.resize(reliability.size());
            I.resize(reliability.size());
            zValues.resize(reliability.size());
        };
        void getReliability(std::vector<double>& _reliability) const { _reliability = reliability; };
    };
}

#endif // _I_INDEX_SCHEME_OPT_METHOD_H_
