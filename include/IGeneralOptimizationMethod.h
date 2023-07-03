#ifndef I_GENERAL_OPTIMIZATION_METHOD_H_
#define I_GENERAL_OPTIMIZATION_METHOD_H_

#include <vector>

using std::vector;

class IGeneralOptimizationMethod {
protected:
    int dimension;
    vector<double> lowerBound, upBound; // area of search

    virtual void solveImplementation() = 0;
    virtual bool solveTestImplementation() = 0;

public:
    IGeneralOptimizationMethod(int _dimension, const vector<double> &_lowerBound, const vector<double> &_upBound)
                               : dimension(_dimension), lowerBound(_lowerBound), upBound(_upBound) {}

    void setDimension(int _dimension) { dimension = _dimension; };
    void setLowerBound(const vector<double> &_lowerBound) { lowerBound = _lowerBound; };
    void setUpBound(const vector<double> &_upBound) { upBound = _upBound; };
    void setLowerUpBounds(const vector<double> &_lowerBound, const vector<double> &_upBound) {
        lowerBound = _lowerBound; upBound = _upBound;
    };

    int getDimension() const { return dimension; };
    void getLowerBound(vector<double>& _lowerBound) const { _lowerBound = lowerBound; };
    void getUpBound(vector<double>& _upBound) const { _upBound = upBound; };
    void getLowerUpBounds(vector<double>& _lowerBound, vector<double>& _upBound) const {
        _lowerBound = lowerBound; _upBound = upBound;
    };
};

#endif // I_GENERAL_OPTIMIZATION_METHOD_H_
