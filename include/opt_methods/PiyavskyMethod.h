#ifndef PIYAVSKY_METHOD_H_
#define PIYAVSKY_METHOD_H_

#include <limits>

#include <base_classes/opt_methods/IConstantEstimationOptMethod.h>
#include <opt_methods/ScanningMethod.h>
#include <my_math.h>

static const double epsilon = 1e-14;

template <typename OptProblemType>
class PiyavskyMethod : public ScanningMethod<OptProblemType>,
                       public opt::IConstantEstimationOptMethod
{
public:
    using GeneralMethod = typename ScanningMethod<OptProblemType>::GeneralMethod;

    using typename ScanningMethod<OptProblemType>::Task;
    using typename ScanningMethod<OptProblemType>::StoppingConditions;
    using typename ScanningMethod<OptProblemType>::Result;

    struct Parameters : public ScanningMethod<OptProblemType>::Parameters {
        double reliability;

        Parameters(double _accuracy = 0.001, double _error = 0.001,
                   int _maxTrials = 1000, int _maxFevals = 1000,
                   double _reliability = 2.0):
            ScanningMethod<OptProblemType>::Parameters(_accuracy, _error, _maxTrials, _maxFevals),
            reliability(_reliability)
        {};
    };

protected:
    void estimatingConstant() override;

    void calcCharacteristic() override;

    double selectNewPoint() override;

public:
    PiyavskyMethod(const OneDimensionalProblem &_problem = OneDimensionalProblem(),
                   const Parameters &parameters = Parameters()):
        ScanningMethod<OptProblemType>(_problem, parameters),
        opt::IConstantEstimationOptMethod(parameters.reliability)
    {};

    void setParameters(const typename GeneralMethod::Parameters &parameters) override {
        ScanningMethod<OptProblemType>::setParameters(parameters);
        reliability = static_cast<const Parameters&>(parameters).reliability;
    };
    void getParameters(typename GeneralMethod::Parameters &parameters) const override {
        ScanningMethod<OptProblemType>::getParameters(parameters);
        static_cast<Parameters&>(parameters).reliability = reliability;
    };
};

template <typename OptProblemType>
void PiyavskyMethod<OptProblemType>::estimatingConstant() {
    static double M = 0.0;
    Trial trialPointT;

    // with optimization(const)
    if (this->trialPoints[this->t].x == this->problem.getSearchArea().upBound) {
        M = std::abs((this->trialPoints[1].z - this->trialPoints[0].z) /
                     (this->trialPoints[1].x - this->trialPoints[0].x));
    } else {
        trialPointT = this->trialPoints[this->t];
        M = std::max({ M, std::abs((trialPointT.z - this->trialPoints[(size_t)this->t - 1].z) / 
                                   (trialPointT.x - this->trialPoints[(size_t)this->t - 1].x)), 
                          std::abs((this->trialPoints[(size_t)this->t + 1].z - trialPointT.z) / 
                                   (this->trialPoints[(size_t)this->t + 1].x - trialPointT.x)) });
    }

    // without optimization(linear)
    // double Mtmp;
    // M = 0.0;
    // int trialPointsSize = (int)this->trialPoints.size();
    // for (int i = 1; i < trialPointsSize; i++) {
    //     Mtmp = std::abs(this->trialPoints[i].z - this->trialPoints[(size_t)i - 1].z) /
    //                    (this->trialPoints[i].x - this->trialPoints[(size_t)i - 1].x);
    //     if (Mtmp > M) M = Mtmp;
    // }

    constantEstimation = (M <= epsilon) ? 1.0 : reliability * M;
}

template <typename OptProblemType>
void PiyavskyMethod<OptProblemType>::calcCharacteristic() {
    double R = -std::numeric_limits<double>::infinity(), Rtmp;
    Trial trialPrev = this->trialPoints[0], trialNext;
    
    int trialPointsSize = (int)this->trialPoints.size();
    for (int i = 1; i < trialPointsSize; i++) {
        trialNext = this->trialPoints[i];
        Rtmp = 0.5 * (constantEstimation * (trialNext.x - trialPrev.x) - (trialNext.z + trialPrev.z));
        if (Rtmp > R) {
            R = Rtmp;
            this->t = i;
        }
        trialPrev = trialNext;
    }
}

template <typename OptProblemType>
double PiyavskyMethod<OptProblemType>::selectNewPoint() {
    estimatingConstant();

    calcCharacteristic();

    return 0.5 * ((this->trialPoints[this->t].x + this->trialPoints[(size_t)this->t - 1].x) -
                  (this->trialPoints[this->t].z - this->trialPoints[(size_t)this->t - 1].z) / constantEstimation);
}

#endif // PIYAVSKY_METHOD_H_
