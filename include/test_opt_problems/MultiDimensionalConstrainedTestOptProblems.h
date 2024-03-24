#ifndef _MULTI_DIMENSIONAL_CONSTRAINED_TEST_OPT_PROBLEMS_H_
#define _MULTI_DIMENSIONAL_CONSTRAINED_TEST_OPT_PROBLEMS_H_

#include <vector>
#include <string>
#include <limits>

#include <opt_problems/ConstrainedOptProblem.h>
#include <general/structures/search_areas/MultiDimensionalSearchArea.h>
#include <MyMath.h>

const int numberBlocks = 2;
const int numberFunctions[2] = { 4, 2 };
const std::vector<std::string> blockNames{ "sample", "test" };

// TODO: replace with sampleProblems
const std::vector<MultiDimensionalConstrainedOptProblem> sampleTasks{
    MultiDimensionalConstrainedOptProblem(
        [] (std::vector<double> x, size_t j) {
            switch (j) {
                case 0: return 1.0 - x[0] - x[1];
                default: return std::numeric_limits<double>::quiet_NaN();
            }
        },
        blockNames[0], 1, 0,
        opt::MultiDimensionalSearchArea(2, std::vector<double>{ -4.0, -4.0 }, std::vector<double>{ 4.0, 4.0 }),
        std::vector<std::vector<double>>{ std::vector<double>{ 4.0, 4.0 } }, -7.0)
        // TODO: add lipshitzConst
};

/* const std::vector<OneDimensionalConstrainedOptProblem> testTasks{
    OneDimensionalConstrainedOptProblem(
        [] (std::vector<double> x, int j) {
            switch (j) {
                case 0: return std::exp(-std::sin(3.0 * x)) - 1.0 / 10.0 * std::pow(x - 1.0 / 2.0, 2) - 1.0;
                case 1: return -13.0 / 6.0 * x + std::sin(13.0 / 4.0 * (2.0 * x + 5.0)) - 53.0 / 12.0;
                default: return std::numeric_limits<double>::quiet_NaN();
            }
        },
        blockNames[1], 1, 1,
        opt::OneDimensionalSearchArea(-2.5, 1.5), std::vector<double>{ 1.05738 }, -7.61284448, 8.666667,
        std::vector<double>{ 4.640837 })
} */

/* double f2Sample(std::vector<double> x, int j) {
    switch (j) {
        case 0: return (x[0] - 1.0) * (x[0] - 1.0) / 5.0 + (x[1] - 1.0) * (x[1] - 1.0) / 5.0;
        default: return std::numeric_limits<double>::quiet_NaN();
    }
}

double f3Sample(std::vector<double> x, int j) {
    switch (j) {
        case 0: return 1.0 - x[0] - x[1];
        case 1: return x[0] * x[0] / 5.0 + x[1] * x[1] / 5.0;
        default: return std::numeric_limits<double>::quiet_NaN();
    }
}

double f4Sample(std::vector<double> x, int j) {
    switch (j) {
        case 0: return (x[0] - 2.0) * (x[0] - 2.0) + (x[1] - 2.0) * (x[1] - 2.0) - 2.0;
        case 1: return x[0] * x[0] / 5.0 + x[1] * x[1] / 5.0;
        default: return std::numeric_limits<double>::quiet_NaN();
    }
} */

/* double f1Test(std::vector<double> x, int j) {
    switch (j) {
        case 0: return 0.01 * (std::pow((x[0] - 2.2), 2) + std::pow((x[1] - 1.2), 2) - 2.25);
        case 1: return 100.0 * (1.0 - std::pow((x[0] - 2.0), 2) / 1.44 - std::pow(0.5 * x[1], 2));
        case 2: return 10.0 * (x[1] - 1.5 - 1.5 * std::sin(6.283 * (x[0] - 1.75)));
        case 3: return -1.5 * x[0] * x[0] * exp(1.0 - x[0] * x[0] - 20.25 * std::pow((x[0] - x[1]), 2)) -
                       std::pow(0.5 * (x[0] - 1.0) * (x[1] - 1.0), 4) * std::exp(2.0 - std::pow(0.5 * (x[0] - 1.0), 4) -
                       std::pow(x[1] - 1.0, 4));
        default: return std::numeric_limits<double>::quiet_NaN();
    }
} */

/*
const double C[20] = { 75.1963666677, -3.8112755343, 0.1269366345, -0.0020567665, 0.000010345,
                      -6.8306567631, 0.0302344793, -0.0012813448, 0.0000352559, -0.0000002266,
                       0.2564581253, -0.0034604030, 0.0000135139, -28.1064434908, -0.0000052375,
                      -0.0000000063, 0.0000000007, 0.0003405462, -0.0000016638, -2.8673112392 };

double f2Test(std::vector<double> x, int j) {
    switch (j) {
        case 0: return 450.0 - x[0] * x[1];
        case 1: return (0.1 * x[0] - 1.0) * (0.1 * x[0] - 1.0) - x[1];
        case 2: return 8.0 * (x[0] - 40.0) - (x[1] - 30.0) * (x[1] - 55.0);
        case 3: return x[1] + (x[0] - 35.0) * (x[0] - 30.0) / 125.0 - 80.0;
        case 4: return -(C[0] + C[1] * x[0] + C[2] * x[0] * x[0] + C[3] * std::pow(x[0], 3) + C[4] * std::pow(x[0], 4) + C[5] * x[1] +
                         C[6] * x[0] * x[1] + C[7] * x[0] * x[0] * x[1] + C[8] * std::pow(x[0], 3) * x[1] + C[9] * std::pow(x[0], 4) * x[1] +
                         C[10] * x[1] * x[1] + C[11] * std::pow(x[1], 3) + C[12] * std::pow(x[1], 4) + C[13] / (x[1] + 1) +
                         C[14] * x[0] * x[0] * x[1] * x[1] + C[15] * std::pow(x[0], 3) * x[1] * x[1] +
                         C[16] * std::pow(x[0], 3) * std::pow(x[1], 3) + C[17] * x[0] * x[1] * x[1] +
                         C[18] * x[0] * std::pow(x[1], 3) + C[19] * std::exp(0.0005 * x[0] * x[1]));
        default: return std::numeric_limits<double>::quiet_NaN();
    }
} */

#endif // _MULTI_DIM_CONSTRAINED_TEST_PROBLEMS_H_
