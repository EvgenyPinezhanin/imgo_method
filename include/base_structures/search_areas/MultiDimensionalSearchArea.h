#ifndef MULTI_DIMENSIONAL_SEARCH_AREA_H_
#define MULTI_DIMENSIONAL_SEARCH_AREA_H_

#include <vector>

using std::vector;

namespace opt {
    class MultiDimensionalSearchArea {
        int dimension;
        vector<double> lowerBound, upBound;

        MultiDimensionalSearchArea(int _dimension = 2, const vector<double> &_lowerBound = vector<double>(),
                                   const vector<double> &_upBound = vector<double>())
                                  : dimension(_dimension), lowerBound(_lowerBound), upBound(_upBound) {};
    };
}

#endif // MULTI_DIMENSIONAL_SEARCH_AREA_H_
