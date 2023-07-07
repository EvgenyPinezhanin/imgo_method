#ifndef ONE_DIMENSIONAL_SEARCH_AREA_H_
#define ONE_DIMENSIONAL_SEARCH_AREA_H_

class OneDimensionalSearchArea {
protected:
    double lowerBound, upBound;

public:
    OneDimensionalSearchArea(double _lowerBound, double _upBound) : lowerBound(_lowerBound), upBound(_upBound) {};

    void setLowerBound(double _lowerBound) { lowerBound = _lowerBound; };
    double getLowerBound() const { return lowerBound; };

    void setUpBound(double _upBound) { upBound = _upBound; };
    double getUpBound() const { return upBound; };

    void setBounds(double _lowerBound, double &_upBound) { lowerBound = _lowerBound; upBound = _upBound; };
    void getBounds(double &_lowerBound, double &_upBound) const { _lowerBound = lowerBound; _upBound = upBound; };
};

#endif // ONE_DIMENSIONAL_SEARCH_AREA_H_
