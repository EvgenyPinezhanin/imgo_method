#ifndef SCANNING_METHOD_H_
#define SCANNING_METHOD_H_

#include <base_classes/opt_methods/ICharacteristicOptMethod.h>
#include <base_classes/opt_methods/IGeneralNumericalOptMethod.h>
#include <opt_problems/OneDimensionalProblem.h>
#include <trials/Trial.h>

class ScanningMethod : public opt::ICharacteristicOptMethod<Trial>,
    public opt::IGeneralNumericalOptMethod<Trial, OneDimensionalProblem>
{
public:
    using GeneralNumMethod = IGeneralNumericalOptMethod<Trial, OneDimensionalProblem>;

    struct Task : public GeneralNumMethod::Task {
        int blockNumber, functionNumber;

        Task(const std::string &_name, int _blockNumber, int _functionNumber,
             const OneDimensionalProblem &_problem, Parameters &_parameters, bool _use = true):
            GeneralNumMethod::Task(_name, _problem, _parameters, _use),
            blockNumber(_blockNumber),
            functionNumber(_functionNumber)
        {};
    };

    class Report : public GeneralNumMethod::IReport {
    protected:
        void printOptProblem(std::ostream &stream, const OneDimensionalProblem &optProblem) const override;
        void printErrorEstimate(std::ostream &stream, const OneDimensionalProblem &optProblem,
                                const GeneralMethod::Result &result) const override;

    public:
        Report():
            GeneralNumMethod::IReport()
        {};
    };
private:
    int cycleBoundary;

protected:
    void insertInSorted(const Trial &trial) override;

    void calcCharacteristic() override;

    Trial newTrial(const OneDimensionalProblem::Point &x) override;
    OneDimensionalProblem::Point selectNewPoint() override;

    double estimateSolution(OneDimensionalProblem::Point &x) const override;

    bool stopConditions() override;
    bool stopConditionsTest() override;

public:
    ScanningMethod(const OneDimensionalProblem &_problem = OneDimensionalProblem(),
                   const Parameters &parameters = Parameters()):
        opt::ICharacteristicOptMethod<Trial>(),
        GeneralNumMethod(_problem, parameters),
        cycleBoundary(0)
    {};

    void setParameters(const GeneralMethod::Parameters &parameters) override {
        GeneralNumMethod::setParameters(parameters);
    };
    void getParameters(GeneralMethod::Parameters &parameters) const override {
        GeneralNumMethod::getParameters(parameters);
    };

    void solve(GeneralMethod::Result &result) override;
    bool solveTest(GeneralMethod::Result &result) override;
};

#endif // SCANNING_METHOD_H_
