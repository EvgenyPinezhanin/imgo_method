#ifndef GENERAL_RESULT_METHOD_H_
#define GENERAL_RESULT_METHOD_H_

namespace om {
    template <typename PointType>
    struct GeneralResultMethod {
        PointType point;
        double value;

        GeneralResultMethod(const PointType _point = PointType(), double _value = 0.0)
                           : point(_point), value(_value) {};
    };
}

#endif // GENERAL_RESULT_METHOD_H_
