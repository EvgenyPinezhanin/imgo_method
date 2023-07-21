#ifndef I_CONSTANT_ESTIMATION_OPT_METHOD_H_
#define I_CONSTANT_ESTIMATION_OPT_METHOD_H_

namespace om {
    class IConstantEstimationOptMethod {
    protected:
        double reliability; 
        double constantEstimation;

    public:
        IConstantEstimationOptMethod(double _reliability) : reliability(_reliability), constantEstimation(0.0) {};

        void setReliability(double _reliability) { reliability = _reliability; };
        double getReliability() const { return reliability; };

        double getConstantEstimation() const { return constantEstimation; };
    };
}

#endif // I_CONSTANT_ESTIMATION_OPT_METHOD_H_
