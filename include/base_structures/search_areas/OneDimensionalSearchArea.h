#ifndef ONE_DIMENSIONAL_SEARCH_AREA_H_
#define ONE_DIMENSIONAL_SEARCH_AREA_H_

namespace opt {
    struct OneDimensionalSearchArea {
        double lowerBound, upBound;

        OneDimensionalSearchArea(double _lowerBound = 0.0, double _upBound = 1.0)
                                : lowerBound(_lowerBound), upBound(_upBound) {};
    };
}

#endif // ONE_DIMENSIONAL_SEARCH_AREA_H_
