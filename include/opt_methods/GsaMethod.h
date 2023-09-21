#ifndef GSA_METHOD_H_
#define GSA_METHOD_H_

#include <opt_methods/PiyavskyMethod.h>

class GsaMethod : public PiyavskyMethod {
public:
    using Parameters = PiyavskyMethod::Parameters;

protected:
    void calcCharacteristic() override;

public:
    GsaMethod(const OneDimensionalProblem &_problem = OneDimensionalProblem(),
              const Parameters &parameters = Parameters()):
        PiyavskyMethod(_problem, parameters)
    {};
};

#endif // GSA_METHOD_H_
