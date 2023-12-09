#ifndef IMGO_METHOD_H_
#define IMGO_METHOD_H_

#include <vector>
#include <algorithm>
#include <limits>
#include <iterator>

#include <base_classes/opt_methods/IGeneralNumericalOptMethod.h>
#include <base_classes/opt_methods/ICharacteristicOptMethod.h>
#include <base_classes/opt_methods/IIndexSchemeOptMethod.h>
#include <base_structures/trials/IndexTrial.h>
#include <opt_problems/OneDimConstrainedProblem.h>
#include <my_math.h>

static const double epsilon = 1e-14;

template <typename OptProblemType>
class ImgoMethod : public opt::IGeneralNumericalOptMethod<opt::IndexTrial, OptProblemType>, 
    public opt::ICharacteristicOptMethod<opt::IndexTrial>, public opt::IIndexSchemeOptMethod
{
public:
    using GeneralNumMethod = opt::IGeneralNumericalOptMethod<opt::IndexTrial, OptProblemType>;
    using typename GeneralNumMethod::GeneralMethod;

    using typename GeneralNumMethod::StoppingConditions;
    using typename GeneralNumMethod::Result;

    struct Parameters : public GeneralNumMethod::Parameters {
        std::vector<double> reliability;
        double d;

        Parameters(double _accuracy = 0.001, double _error = 0.001,
                   int _maxTrials = 1000, int _maxFevals = 1000,
                   std::vector<double> _reliability = 2.0, double _d = 0.0):
            GeneralNumMethod::Parameters(_accuracy, _error, _maxTrials, _maxFevals),
            reliability(_reliability), d(_d)
        {};
    };

    struct Task : public GeneralNumMethod::Task {
        std::string blockName;
        int functionNumber;

        Task(const std::string &_name, std::string _blockName, int _functionNumber,
             const OptProblemType &_problem, Parameters &_parameters, bool _use = true):
            GeneralNumMethod::Task(_name, _problem, _parameters, _use),
            blockName(_blockName),
            functionNumber(_functionNumber)
        {};
    };

    class Report : public GeneralNumMethod::IReport {
    protected:
        void printOptProblem(std::ostream &stream, const OptProblemType &optProblem) const override;
        void printMethodParameters(std::ostream &stream,
                                   const typename GeneralMethod::Parameters &parameters) const override;
        void printResultMethod(std::ostream &stream,
                               const typename GeneralMethod::Result &result) const override;
        void printErrorEstimate(std::ostream &stream, const OptProblemType &optProblem,
                                const typename GeneralMethod::Result &result) const override;

    public:
        Report():
            GeneralNumMethod::IReport()
        {};
    };

protected:
    void insertInSorted(const Trial &trial) override;

    void calcCharacteristic() override;

    Trial newTrial(const typename OptProblemType::Point &x) override;
    typename OptProblemType::Point selectNewPoint() override;

    double estimateSolution(typename OptProblemType::Point &x) const override;

    bool stopConditions() override;
    bool stopConditionsTest() override;

    virtual void estimatingConstants() = 0;

    virtual void calcI() = 0;
    virtual void calcZValues() = 0;

    using GeneralNumMethod::setResult;
private:
    function<double(double, int)> f;
    double r, d;

    TrialConstrained lastTrial;
    int lastTrialPos;

    vector<bool> IIsCalc;
    vector<double> mu;
    vector<double> zStar;

    TrialConstrained newTrial(double x) override;
    double newPoint(int t) override;
    double selectNewPoint(int &t) override;

public:
    ImgoMethod(function<double(double, int)> _f = nullptr, int _numberConstraints = 0, double _a = 0.0, double _b = 1.0, 
               double _r = 2.0, double _d = 0.0, double _eps = 0.0001, int _maxTrials = 1000, int _maxFevals = 1000)
               : OptimizationMethodConstrained(nullptr, 1, _numberConstraints, vector<double>{_a}, vector<double>{_b},
               _eps, _maxTrials, _maxFevals), f(_f), r(_r), d(_d), lastTrial(0.0, 0.0, 0), lastTrialPos(0),
               I((size_t)numberConstraints + 1), calcI((size_t)numberConstraints + 1), mu((size_t)numberConstraints + 1),
               zStar((size_t)numberConstraints + 1) {}

    void setParameters(const typename GeneralMethod::Parameters &parameters) override {
        GeneralNumMethod::setParameters(parameters);
    };
    void getParameters(typename GeneralMethod::Parameters &parameters) const override {
        GeneralNumMethod::getParameters(parameters);
    };

    void solve(typename GeneralMethod::Result &result) override;
    bool solveTest(typename GeneralMethod::Result &result) override;

    void setF(const function<double(double, int)> &_f) { f = _f; };
    void setA(double _a) { OptimizationMethod::setA(vector<double>{ _a }); };
    void setB(double _b) { OptimizationMethod::setB(vector<double>{ _b }); };
    void setAB(double _a, double _b) { OptimizationMethod::setAB(vector<double>{ _a }, vector<double>{ _b }); };
    void setNumberConstraints(int _numberConstraints);
    void setR(double _r) { r = _r; };
    void setD(double _d) { d = _d; };

    function<double(double, int)> getF() const { return f; };
    double getA() const { return A[0]; };
    double getB() const { return B[0]; };
    double getR() const { return r; };
    double getD() const { return d; };

    void getL(vector<double> &L) const { L = mu; };

    bool solveTest(double xOpt, int &numberTrials, int &numberFevals);

        ScanningMethod(const OptProblemType &_problem = OptProblemType(),
                   const Parameters &parameters = Parameters()):
        opt::ICharacteristicOptMethod<Trial>(),
        GeneralNumMethod(_problem, parameters),
        cycleBoundary(0)
    {};
};

template <typename OptProblemType>
void ImgoMethod<OptProblemType>::Report::printOptProblem(
    std::ostream &stream,
    const OptProblemType &optProblem) const
{
    stream << "Number of constraints = " << optProblem.numberConstraints << "\n";
    stream << "[a; b] = [" << optProblem.getSearchArea().lowerBound << "; " <<
                              optProblem.getSearchArea().upBound << "]"<< "\n";
    std::vector<double> optimalPoints;
    optProblem.getOptimalPoints(optimalPoints);
    if (!optimalPoints.empty()) {
        stream << "X* = (" << optimalPoints[0];
        int numberOptimalPoints = optimalPoints.size();
        for (int i = 1; i < numberOptimalPoints; i++) {
            stream << "; " << optimalPoints[i];
        }
        stream << ")\n";
        stream << "f(X*) = " << optProblem.getOptimalValue() << "\n";
    }

    // cout << "Lipschitz constant:" << "\n";
    // cout << "L*(f) = " << lipschitzConst[numberConstraints] << "\n";
    // for (int j = 0; j < numberConstraints; j++) {
    //     cout << "L*(g" << j + 1 << ") = " << lipschitzConst[j] << "\n";
    // }
}

template <typename OptProblemType>
void ImgoMethod<OptProblemType>::Report::printMethodParameters(
    std::ostream &stream,
    const typename GeneralMethod::Parameters &parameters) const
{
    GeneralNumMethod::Report::printMethodParameters(stream, parameters);
    auto parametersCast = static_cast<const Parameters&>(parameters);

    stream << "reliability = (" << parametersCast.reliability << "\n";
    size_t numberConstraints = parametersCast.reliability.size();
    for (size_t i = 0; i < numberConstraints - 1; ++i) {
        stream << parametersCast.reliability[i] << ", ";
    }
    stream << parametersCast.reliability[numberConstraints - 1] << ")" << "\n";
    stream << "d = " << parametersCast.accuracy << "\n";
}

template <typename OptProblemType>
void ImgoMethod<OptProblemType>::Report::printResultMethod(
    std::ostream &stream,
    const typename GeneralMethod::Result &result) const
{
    GeneralNumMethod::Report::printResultMethod(stream, result);
    
    // cout << "Estimation of the Lipschitz constant:" << "\n";
    // cout << "L(f) = " << estLipschitzConst[numberConstraints] << "\n";
    // for (int j = 0; j < numberConstraints; j++) {
    //     cout << "L(g" << j + 1 << ") = " << estLipschitzConst[j] << "\n";
    // }
}

template <typename OptProblemType>
void ImgoMethod<OptProblemType>::Report::printErrorEstimate(
    std::ostream &stream,
    const OptProblemType &optProblem,
    const typename GeneralMethod::Result &result) const
{
    double point = result.point;
    std::vector<double> optimalPoints;
    optProblem.getOptimalPoints(optimalPoints);
    if (!optimalPoints.empty()) {
        auto iter = std::min_element(optimalPoints.begin(), optimalPoints.end(),
            [&point] (const double &point1, const double &point2) {
                return std::abs(point1 - point) < std::abs(point2 - point);
            });
        stream << "|X* - X| = " << std::abs(*iter - point) << "\n";
        stream << "|f(X*) - f(X)| = " << std::abs(optProblem.getOptimalValue() - result.value) << "\n";
    }
}

template <typename OptProblemType>
void ScanningMethod<OptProblemType>::solve(typename GeneralMethod::Result &result) {
    this->trialPoints.clear();
    this->numberFevals = 0;

    this->trialPoints.push_back(newTrial(this->problem.getSearchArea().lowerBound));
    this->trialPoints.push_back(newTrial(this->problem.getSearchArea().upBound));
    t = 1;

    this->numberTrials = 2;
    cycleBoundary = 0;

    Trial trial;
    double xNew;
    while(!stopConditions()) {
        xNew = selectNewPoint();

        trial = newTrial(xNew);
        this->numberTrials++;

        insertInSorted(trial);
    }

    setResult(result);
}

template <typename OptProblemType>
void ImgoMethod<OptProblemType>::solve(typename GeneralMethod::Result &result) {
    for (size_t i = 0; i < this->problem.numberConstraints; ++i) {
        this->I[i].clear();
        IIsCalc[i] = false;
    }
    this->trialPoints.clear();
    this->numberFevals = 0;

    this->trialPoints.push_back(newTrial(this->problem.getSearchArea().lowerBound));
    this->trialPoints.push_back(newTrial(this->problem.getSearchArea().upBound));
    t = 1;

    lastTrial = newTrial(A[0]);
    trialPoints.push_back(lastTrial);
    insertInSorted(I[(size_t)lastTrial.nu - 1], lastTrial);
    lastTrial = newTrial(B[0]);
    trialPoints.push_back(lastTrial);
    lastTrialPos = insertInSorted(I[(size_t)lastTrial.nu - 1], lastTrial);
    numberTrials = 2;

    double xNew;
    int t;
    while(true) {
        numberTrials++;

        // Steps 3, 4, 5, 6, 7
        xNew = selectNewPoint(t);
        lastTrial = newTrial(xNew);

        // Step 1
        insertInSorted(trialPoints, lastTrial);

        // Step 2
        lastTrialPos = insertInSorted(I[(size_t)lastTrial.nu - 1], lastTrial);

        // Stop conditions
        if (trialPoints[t].x - trialPoints[(size_t)t - 1].x <= eps) break;
        if (this->numberFevals >= maxFevals || numberTrials >= maxTrials) break;
    }
    numberFevals = this->numberFevals;
    x = searchMin(trialPoints, numberConstraints);
}

bool ImgoMethod::solveTest(double xOpt, int &numberTrials, int &numberFevals) {
    for (int nu = 0; nu < numberConstraints + 1; nu++) {
        I[nu].clear();
        calcI[nu] = false;
    }
    trialPoints.clear();
    this->numberFevals = 0;

    lastTrial = newTrial(A[0]);
    trialPoints.push_back(lastTrial);
    insertInSorted(I[(size_t)lastTrial.nu - 1], lastTrial);
    lastTrial = newTrial(B[0]);
    trialPoints.push_back(lastTrial);
    lastTrialPos = insertInSorted(I[(size_t)lastTrial.nu - 1], lastTrial);
    numberTrials = 2;

    double xNew;
    int t;
    while (true) {
        numberTrials++;

        // Steps 3, 4, 5, 6, 7
        xNew = selectNewPoint(t);
        lastTrial = newTrial(xNew);

        // Step 1
        insertInSorted(trialPoints, lastTrial);

        // Step 2
        lastTrialPos = insertInSorted(I[(size_t)lastTrial.nu - 1], lastTrial);

        // Stop conditions
        if (abs(xNew - xOpt) <= eps) {
            numberFevals = this->numberFevals;
            return true;
        }
        if (this->numberFevals >= maxFevals || numberTrials >= maxTrials) {
            numberFevals = this->numberFevals;
            return false;
        }
    }
}


#endif // IMGO_METHOD_H_

inline int insertInSorted(vector<TrialConstrained> &trials, TrialConstrained trial) {
    vector<TrialConstrained>::iterator iter = trials.begin();
    vector<TrialConstrained>::iterator iterEnd = trials.end();
    int pos = 0;
    while(true) {
        if (iter == iterEnd || iter->x > trial.x) break;
        iter++; pos++;
    }
    trials.insert(iter, trial);
    return pos;

        vector<Trial>::iterator iter = trialPoints.begin();
    int dist = trialPoints.size();

    advance(iter, dist / 2);
    while (true) {
        dist /= 2;
        if (trial.x < iter->x) advance(iter, -dist / 2);
        else {
            if (trial.x < (iter + 1)->x) break;
            else advance(iter, dist / 2);
        }
    }
    iter = trialPoints.insert(iter, trial);

    return distance(trialPoints.begin(), iter);
}

inline double searchMin(vector<TrialConstrained> &trials, int numberConstraints) {
    double z = numeric_limits<double>::infinity(), x = 0.0;
    size_t sizeTrials = trials.size();
    for (int i = 0; i < sizeTrials; i++) {
        if (trials[i].nu == numberConstraints + 1 && trials[i].z < z) {
            z = trials[i].z;
            x = trials[i].x;
        }
    }
    return x;
}

TrialConstrained ImgoMethod::newTrial(double x) {
    TrialConstrained trial(x);
    for (int j = 1; j <= numberConstraints + 1; j++) {
        numberFevals++;
        if ((f(x, j) > 0) || (j == numberConstraints + 1)) {
            trial.z = f(x, j);
            trial.nu = j;
            break;
        }
    }
    return trial;
}

double ImgoMethod::newPoint(int t) {
    if (trialPoints[t].nu != trialPoints[(size_t)t - 1].nu) {
        return (trialPoints[t].x + trialPoints[(size_t)t - 1].x) / 2.0;
    } else {
        return (trialPoints[t].x + trialPoints[(size_t)t - 1].x) / 2.0 - 
               (trialPoints[t].z - trialPoints[(size_t)t - 1].z) / (2.0 * r * mu[(size_t)trialPoints[t].nu - 1]);
    }
}

double ImgoMethod::selectNewPoint(int &t) {
    int nuLastTrial;
    size_t sizeI;
    double muTmp;

    // Step 3
    // with optimization(const)
    nuLastTrial = lastTrial.nu - 1;
    sizeI = I[nuLastTrial].size();
    for (int nu = 0; nu < numberConstraints + 1; nu++) {
        if (!calcI[nu]) mu[nu] = 0.0;
    }
    if (I[nuLastTrial].size() >= 3) {
        if (lastTrialPos == 0) {
            mu[nuLastTrial] = max({ mu[nuLastTrial], abs(I[nuLastTrial][1].z - I[nuLastTrial][0].z) /
                                                     pow(I[nuLastTrial][1].x - I[nuLastTrial][0].x, 1.0 / n) });  
        } else if (lastTrialPos == I[nuLastTrial].size() - 1) {
            mu[nuLastTrial] = max({ mu[nuLastTrial],
                                    abs(I[nuLastTrial][sizeI - 1].z - I[nuLastTrial][sizeI - 2].z) / 
                                    pow(I[nuLastTrial][sizeI - 1].x - I[nuLastTrial][sizeI - 2].x, 1.0 / n) });
        } else {
            mu[nuLastTrial] = max({ mu[nuLastTrial],
                                    abs(I[nuLastTrial][lastTrialPos].z - I[nuLastTrial][(size_t)lastTrialPos - 1].z) / 
                                    pow(I[nuLastTrial][lastTrialPos].x - I[nuLastTrial][(size_t)lastTrialPos - 1].x, 1.0 / n),
                                    abs(I[nuLastTrial][(size_t)lastTrialPos + 1].z - I[nuLastTrial][lastTrialPos].z) / 
                                    pow(I[nuLastTrial][(size_t)lastTrialPos + 1].x - I[nuLastTrial][lastTrialPos].x, 1.0 / n) });
        }
    } else if (I[nuLastTrial].size() == 2) {
        mu[nuLastTrial] = max({ mu[nuLastTrial], abs(I[nuLastTrial][1].z - I[nuLastTrial][0].z) /
                                                 pow(I[nuLastTrial][1].x - I[nuLastTrial][0].x, 1.0 / n) });
    }
    if (abs(mu[nuLastTrial]) > epsilon) calcI[nuLastTrial] = true;
    for (int nu = 0; nu < numberConstraints + 1; nu++) {
        if (abs(mu[nu]) <= epsilon) mu[nu] = 1.0;
    }

    // with optimization(linear)
    // nuLastTrial = lastTrial.nu - 1;
    // sizeI = I[nuLastTrial].size();
    // for (int nu = 0; nu < numberConstraints + 1; nu++) {
    //     if (!calcI[nuLastTrial]) mu[nu] = 0.0;
    // }
    // for (int i = 1; i < sizeI; i++) {
    //     muTmp = abs(I[nuLastTrial][i].z - I[nuLastTrial][(size_t)i - 1].z) /
    //                (I[nuLastTrial][i].x - I[nuLastTrial][(size_t)i - 1].x);
    //     if (muTmp > mu[nuLastTrial]) {
    //         mu[nuLastTrial] = muTmp;
    //         if (abs(mu[nuLastTrial]) > epsilon) calcI[nuLastTrial] = true;
    //     }
    // }
    // for (int nu = 0; nu < numberConstraints + 1; nu++) {
    //     if (abs(mu[nu]) <= epsilon) mu[nu] = 1.0;
    // }

    // without optimization
    // for (int nu = 0; nu < numberConstraints + 1; nu++) {
    //     mu[nu] = 0.0;
    // }
    // for (int nu = 0; nu < numberConstraints + 1; nu++) {
    //     sizeI = I[nu].size();
    //     for (int i = 1; i < sizeI; i++) {
    //         for (int j = 0; j < i; j++) {
    //             muTmp = abs(I[nu][i].z - I[nu][j].z) / (I[nu][i].x - I[nu][j].x);
    //             if (muTmp > mu[nu]) {
    //                 mu[nu] = muTmp;
    //             }
    //         }
    //     }
    //     if (abs(mu[nu]) <= epsilon) {
    //         mu[nu] = 1.0;
    //     };
    // }

    // Step 4
    for (int nu = 0; nu < numberConstraints + 1; nu++) {
        if (I[nu].size() != 0) {
            zStar[nu] = I[nu][0].z;
            sizeI = I[nu].size();
            for (int i = 1; i < sizeI; i++) {
                if (I[nu][i].z < zStar[nu]) {
                    zStar[nu] = I[nu][i].z;
                }
            }
            if (zStar[nu] <= 0.0 && nu != numberConstraints) {
                zStar[nu] = -mu[nu] * d;
            }
        }
    }

    // Steps 5, 6
    double R = -numeric_limits<double>::infinity(), Rtmp;
    double muV, zStarV, dx;

    size_t sizeTrialPoints = trialPoints.size();
    for (size_t i = 1; i < sizeTrialPoints; i++) {
        dx = trialPoints[i].x - trialPoints[i - 1].x;
        if (trialPoints[i].nu == trialPoints[i - 1].nu) {
            muV = mu[(size_t)trialPoints[i].nu - 1];
            zStarV = zStar[(size_t)trialPoints[i].nu - 1];
            Rtmp = dx + pow(trialPoints[i].z - trialPoints[i - 1].z, 2) / (r * r * muV * muV * dx) -
                   2.0 * (trialPoints[i].z + trialPoints[i - 1].z - 2.0 * zStarV) / (r * muV);
        } else if (trialPoints[i - 1].nu < trialPoints[i].nu) {
            muV = mu[(size_t)trialPoints[i].nu - 1];
            zStarV = zStar[(size_t)trialPoints[i].nu - 1];
            Rtmp = 2.0 * dx - 4.0 * (trialPoints[i].z - zStarV) / (r * muV);
        } else {
            muV = mu[(size_t)trialPoints[i - 1].nu - 1];
            zStarV = zStar[(size_t)trialPoints[i - 1].nu - 1];
            Rtmp = 2.0 * dx - 4.0 * (trialPoints[i - 1].z - zStarV) / (r * muV);
        }
        if (Rtmp > R) {
            R = Rtmp;
            t = (int)i;
        }
    }

    // Step 7
    return newPoint(t);
}

void ImgoMethod::setNumberConstraints(int _numberConstraints) {
    OptimizationMethodConstrained::setNumberConstraints(_numberConstraints);
    I.resize((size_t)numberConstraints + 1);
    calcI.resize((size_t)numberConstraints + 1);
    mu.resize((size_t)numberConstraints + 1);
    zStar.resize((size_t)numberConstraints + 1);
}
