#ifndef MULTI_DIMENSIONAL_SEARCH_AREA_H_
#define MULTI_DIMENSIONAL_SEARCH_AREA_H_

#include <vector>

using std::vector;

class MultiDimensionalSearchArea {
protected:
    int dimension;
    vector<double> lowerBound, upBound;

public:
    MultiDimensionalSearchArea(int _dimension, const vector<double> &_lowerBound, const vector<double> &_upBound)
                               : dimension(_dimension), lowerBound(_lowerBound), upBound(_upBound) {};

    void setDimension(int _dimension) { _dimension = dimension; };
    int getDimension() { return dimension; };

    void setLowerBound(const vector<double> &_lowerBound) { lowerBound = _lowerBound; };
    void getLowerBound(vector<double>& _lowerBound) const { _lowerBound = lowerBound; };

    void setUpBound(const vector<double> &_upBound) { upBound = _upBound; };
    void getUpBound(vector<double>& _upBound) const { _upBound = upBound; };

    void setBounds(const vector<double> &_lowerBound, const vector<double> &_upBound) {
        lowerBound = _lowerBound; upBound = _upBound;
    };
    void getBounds(vector<double>& _lowerBound, vector<double>& _upBound) const {
        _lowerBound = lowerBound; _upBound = upBound;
    };
};

#endif // MULTI_DIMENSIONAL_SEARCH_AREA_H_
