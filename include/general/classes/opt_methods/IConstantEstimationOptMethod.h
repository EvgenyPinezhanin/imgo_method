#ifndef _I_CONSTANT_ESTIMATION_OPT_METHOD_H_
#define _I_CONSTANT_ESTIMATION_OPT_METHOD_H_

namespace opt {
    class IConstantEstimationOptMethod {
    protected:
        double reliability;
        double constantEstimation;

        virtual void estimatingConstant() = 0;

    public:
        IConstantEstimationOptMethod(double _reliability) : reliability(_reliability), constantEstimation(0.0) {};

        void setReliability(double _reliability) { reliability = _reliability; };
        double getReliability() const { return reliability; };
    };
}

#endif // _I_CONSTANT_ESTIMATION_OPT_METHOD_H_
