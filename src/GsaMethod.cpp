#include <opt_methods/GsaMethod.h>

#include <limits>

#include <my_math.h>

using std::numeric_limits;

void GsaMethod::calcCharacteristic() {
    double R = -numeric_limits<double>::infinity(), Rtmp;
    
    int trialPointsSize = (int)trialPoints.size();
    for (int i = 1; i < trialPointsSize; i++) {
        Rtmp = constantEstimation * (trialPoints[i].x - trialPoints[i - 1].x) + (trialPoints[i].z - trialPoints[i - 1].z) *
               (trialPoints[i].z - trialPoints[i - 1].z) / (constantEstimation * (trialPoints[i].x - trialPoints[i - 1].x)) -
               2.0 * (trialPoints[i].z + trialPoints[i - 1].z);
        if (Rtmp > R) {
            R = Rtmp;
            t = i;
        }
    }
}
