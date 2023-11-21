#ifndef GSA_METHOD_H_
#define GSA_METHOD_H_

#include <limits>

#include <opt_methods/PiyavskyMethod.h>
#include <my_math.h>

template <typename OptProblemType>
class GsaMethod : public PiyavskyMethod<OptProblemType> {
public:
    using Parameters = typename PiyavskyMethod<OptProblemType>::Parameters;

protected:
    void calcCharacteristic() override;

public:
    GsaMethod(const OneDimensionalProblem &_problem = OneDimensionalProblem(),
              const Parameters &parameters = Parameters()):
        PiyavskyMethod<OptProblemType>(_problem, parameters)
    {};
};

template <typename OptProblemType>
void GsaMethod<OptProblemType>::calcCharacteristic() {
    double R = -std::numeric_limits<double>::infinity(), Rtmp;
    Trial trialPrev = this->trialPoints[0], trialNext;

    int trialPointsSize = (int)this->trialPoints.size();
    for (int i = 1; i < trialPointsSize; i++) {
        trialNext = this->trialPoints[i];
        Rtmp = this->constantEstimation * (trialNext.x - trialPrev.x) + (trialNext.z - trialPrev.z) *
               (trialNext.z - trialPrev.z) / (this->constantEstimation * (trialNext.x - trialPrev.x)) -
               2.0 * (trialNext.z + trialPrev.z);
        if (Rtmp > R) {
            R = Rtmp;
            this->t = i;
        }
        trialPrev = trialNext;
    }
}

#endif // GSA_METHOD_H_
