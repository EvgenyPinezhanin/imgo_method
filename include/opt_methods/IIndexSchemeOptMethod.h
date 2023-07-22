#ifndef I_INDEX_SCHEME_OPT_METHOD_H_
#define I_INDEX_SCHEME_OPT_METHOD_H_

#include <vector>

#include <trials/IndexTrial.h>

using std::vector;

class IIndexSchemeOptMethod {
protected:
    vector<double> reliability;
    vector<double> constantsEstimation;

/*     vector<vector<IndexTrial>> I;
    vector<double> mu;
    vector<double> zStar; */

public:
    IIndexSchemeOptMethod(const vector<double> &_reliability)
                         : reliability(_reliability), constantsEstimation(0) {};

/*     void setR(const vector<double> &_r) { r = _r; };
    void getR(vector<double>& _r) const { _r = r; };

    void getConstantEstimation(const vector<double> &_constantsEstimation) const { _constantsEstimation = constantsEstimation; }; */
};

#endif // I_INDEX_SCHEME_OPT_METHOD_H_
