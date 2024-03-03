#ifndef ONE_DIMENSIONAL_CONSTRAINED_TEST_OPT_PROBLEMS_H_
#define ONE_DIMENSIONAL_CONSTRAINED_TEST_OPT_PROBLEMS_H_

#include <vector>
#include <string>
#include <limits>

#include <opt_problems/ConstrainedOptProblem.h>
#include <general/structures/search_areas/OneDimensionalSearchArea.h>
#include <MyMath.h>

const int numberBlocks = 2;
const int numberProblems[2] = { 3, 10 };
const std::vector<std::string> blockNames = { "sample", "test" };

const std::vector<OneDimensionalConstrainedOptProblem> sampleTasks{
    OneDimensionalConstrainedOptProblem(
        [] (double x, int j) {
            switch (j) {
                case 0: return std::sin(x);
                default: return std::numeric_limits<double>::quiet_NaN();
            }
        },
        blockNames[0], 1, 0,
        opt::OneDimensionalSearchArea(-4.0, 4.0), std::vector<double>{ -M_PI / 2.0 }, -1.0, 1.0),

    OneDimensionalConstrainedOptProblem(
        [] (double x, int j) {
            switch (j) {
                case 0: return std::sin(x);
                case 1: return -2.0 * x + 3.0;
                default: return std::numeric_limits<double>::quiet_NaN();
            }
        },
        blockNames[0], 2, 1,
        opt::OneDimensionalSearchArea(2.0, 8.0), std::vector<double>{ 2.0 * M_PI }, -9.566371, 2.0,
        std::vector<double>{ 1.0 }),

    OneDimensionalConstrainedOptProblem(
        [] (double x, int j) {
            switch (j) {
                case 0: return x * x - 0.05;
                case 1: return -x + 0.1;
                case 2: return 5.0 * x * x + 3.0 * x - 1.0;
                default: return std::numeric_limits<double>::quiet_NaN();
            }
        },
        blockNames[0], 3, 2,
        opt::OneDimensionalSearchArea(-2.0, 2.0), std::vector<double>{ 0.1 }, -0.65, 23.0,
        std::vector<double>{ 4.0, 1.0 })
};

const std::vector<OneDimensionalConstrainedOptProblem> testTasks{
    OneDimensionalConstrainedOptProblem(
        [] (double x, int j) {
            switch (j) {
                case 0: return std::exp(-std::sin(3.0 * x)) - 1.0 / 10.0 * std::pow(x - 1.0 / 2.0, 2) - 1.0;
                case 1: return -13.0 / 6.0 * x + std::sin(13.0 / 4.0 * (2.0 * x + 5.0)) - 53.0 / 12.0;
                default: return std::numeric_limits<double>::quiet_NaN();
            }
        },
        blockNames[1], 1, 1,
        opt::OneDimensionalSearchArea(-2.5, 1.5), std::vector<double>{ 1.05738 }, -7.61284448, 8.666667,
        std::vector<double>{ 4.640837 }),

    OneDimensionalConstrainedOptProblem(
        [] (double x, int j) {
            switch (j) {
                case 0: return 1.0 / 20.0 - std::exp(-2.0 / 5.0 * (x + 5.0)) * std::sin(4.0 / 5.0 * M_PI * (x + 5.0));
                case 1: return (11.0 * x * x - 10.0 * x + 21.0) / (2.0 * (x * x + 1.0));
                default: return std::numeric_limits<double>::quiet_NaN();
            }
        },
        blockNames[1], 2, 1,
        opt::OneDimensionalSearchArea(-5.0, 5.0), std::vector<double>{ 1.016 }, 5.46063488, 6.372595,
        std::vector<double>{ 2.513269 }),

    OneDimensionalConstrainedOptProblem(
        [] (double x, int j) {
            double sum = 0.0;
            switch (j) {
                case 0: return 3.0 / 2.0 * (std::cos(7.0 / 20.0 * (x + 10.0)) -
                               std::sin(7.0 / 4.0 * (x + 10.0)) + 1.0 / 2.0);
                case 1:
                    for (int i = 1; i <= 5; i++) {
                        sum += std::cos(i * x);
                    }
                    return -sum;
                default: return std::numeric_limits<double>::quiet_NaN();
            }
        },
        blockNames[1], 3, 1,
        opt::OneDimensionalSearchArea(-10.0, 10.0), std::vector<double>{ -5.9921 }, -2.94600839, 13.201241,
        std::vector<double>{ 3.124504 }),

    OneDimensionalConstrainedOptProblem(
        [] (double x, int j) {
            double sum = 0.0;
            switch (j) {
                case 0:
                    for (int i = 1; i <= 5; ++i) {
                        sum += std::cos(5.0 / 4.0 * (i + 1.0) * x + i);
                    }
                    return 6.0 / 25.0 - sum;
                case 1: return 9.0 / 50.0 - 9.0 / 2.0 * std::exp(-(x - 1.0 / 10.0)) * 
                               std::sin(2.0 * M_PI * (x - 1.0 / 10.0));
                case 2: return 4.0 * std::sin(M_PI / 4.0 * x + 1.0 / 20.0) *
                               std::pow(std::pow(std::sin(M_PI / 2.0 * x + 1.0 / 10.0), 3.0) +
                               std::pow(std::cos(M_PI / 2.0 * x + 1.0 / 10.0), 3.0), 2.0);
                default: return std::numeric_limits<double>::quiet_NaN();
            }
        },
        blockNames[1], 4, 2,
        opt::OneDimensionalSearchArea(0.0, 4.0), std::vector<double>{ 2.45956 }, 1.8408089, 12.893183,
        std::vector<double>{ 29.731102, 35.390605 }),

    OneDimensionalConstrainedOptProblem(
        [] (double x, int j) {
            double sum = 0.0;
            switch (j) {
                case 0: return 17.0 / 25.0 - 2.0 / 29763.233 * (-1.0 / 6.0 * std::pow(x, 6) +
                               52.0 / 25.0 * std::pow(x, 5) - 39.0 / 80.0 * std::pow(x, 4) -
                               71.0 / 10.0 * x * x * x + 79.0 / 20.0 * x * x + x - 1.0 / 10.0);
                case 1: return -14.0 / 125.0 * (3.0 * x - 8.0) * std::sin(252.0 / 125.0 * (x + 3.0 / 2.0)) - 1.0 / 2.0;
                case 2: return std::sin(0.423531 * x + 3.13531) + std::sin(10.0 / 3.0 * (0.423531 * x + 3.13531)) +
                               std::log(0.423531 * x + 3.13531) + 0.36634 - 0.355766 * x;
                default: return std::numeric_limits<double>::quiet_NaN();
            }
        },
        blockNames[1], 5, 2,
        opt::OneDimensionalSearchArea(-1.5, 11.0), std::vector<double>{ 9.28491 }, -1.27299809, 2.021595,
        std::vector<double>{ 5.654617, 0.931981 }),

    OneDimensionalConstrainedOptProblem(
        [] (double x, int j) {
            double sum = 0.0;
            switch (j) {
                case 0: return 40.0 * (std::cos(4.0 * x) * (x - std::sin(x)) * std::exp( -(x * x) / 2.0));
                case 1: return 2.0 / 25.0 * (x + 4.0) - std::sin(12.0 / 5.0 * (x + 4.0));
                case 2: return -7.0 / 40.0 * (3.0 * x + 4.0) * std::sin(63.0 / 20.0 * (x + 4.0));
                default: return std::numeric_limits<double>::quiet_NaN();
            }
        },
        blockNames[1], 6, 2,
        opt::OneDimensionalSearchArea(-4.0, 4.0), std::vector<double>{ 2.32396 }, -1.6851399, 8.835339,
        std::vector<double>{ 2.48, 25.108154 }),

    OneDimensionalConstrainedOptProblem(
        [] (double x, int j) {
            double sum = 0.0;
            switch (j) {
                case 0: return std::pow(std::sin(x), 3) * std::exp(-std::sin(3.0 * x)) + 1.0 / 2.0;
                case 1: return std::cos(7.0 / 5.0 * (x + 3.0)) - std::sin(7.0 * (x + 3.0)) + 3.0 / 10.0;
                case 2: return std::exp(-std::cos(4.0 * x - 3.0)) + 1.0 / 250.0 * (4.0 * x - 3.0) * (4.0 * x - 3.0) - 1.0;
                default: return std::numeric_limits<double>::quiet_NaN();
            }
        },
        blockNames[1], 7, 2,
        opt::OneDimensionalSearchArea(-3.0, 2.0), std::vector<double>{ -0.774575 }, -0.47704013, 6.387862,
        std::vector<double>{ 8.332010, 5.359309 }),

    OneDimensionalConstrainedOptProblem(
        [] (double x, int j) {
            double sum = 0.0;
            switch (j) {
                case 0: return std::exp(-std::sin(4.0 * x)) - 1.0 / 10.0 * std::pow(x - 1.0 / 2.0, 2) - 1.0;
                case 1: 
                    for (int i = 1; i <= 5; i++) {
                        sum += std::cos(5.0 * (i + 1.0) * (x + 1.0 / 2.0));
                    }
                    return 3.0 / 10.0 - sum;
                case 2: return (-21.0 / 20.0 * x - 13.0 / 8.0) * std::sin(63.0 / 10.0 * x + 63.0 / 4.0) + 1.0 / 5.0;
                case 3: return std::cos(7.0 / 4.0 * x + 241.0 / 40.0) - std::sin(35.0 / 4.0 * x + 241.0 / 8.0) - 5.0;
                default: return std::numeric_limits<double>::quiet_NaN();
            }
        },
        blockNames[1], 8, 3,
        opt::OneDimensionalSearchArea(-2.5, 1.5), std::vector<double>{ -1.12724 }, -6.60059665, 10.415012,
        std::vector<double>{ 20.184982, 90.598898, 6.372137 }),

    OneDimensionalConstrainedOptProblem(
        [] (double x, int j) {
            double sum = 0.0;
            switch (j) {
                case 0: return 1.0 / 40.0 * (x - 4.0) * (x - 32.0 / 5.0) * (x - 9.0) * (x - 11.0) *
                               std::exp(-1.0 / 10.0 * std::pow(x - 13.0 / 2.0, 2));
                case 1: return (std::pow(std::sin(x + 1.0), 3) + std::pow(std::cos(x + 1.0), 3)) *
                               std::exp(-(x + 1.0) / 10.0);
                case 2: return std::exp(-std::cos(3.0 / 5.0 * (x - 5.0 / 2.0))) +
                               1.0 / 10.0 * std::pow(3.0 / 25.0 * x - 4.0 / 5.0, 2) - 1.0;
                case 3:
                    for (int i = 1; i <= 5; i++) {
                        sum += 1.0 / 5.0 * std::sin((i + 1.0) * x - 1.0) + 2.0;
                    }
                    return sum;
                default: return std::numeric_limits<double>::quiet_NaN();
            }
        },
        blockNames[1], 9, 3,
        opt::OneDimensionalSearchArea(0.0, 14.0), std::vector<double>{ 4.0 }, 9.92218867, 3.843648,
        std::vector<double>{ 0.873861, 1.682731, 1.254588 }),

    OneDimensionalConstrainedOptProblem(
        [] (double x, int j) {
            double sum = 0.0, a, b;
            switch (j) {
                case 0: return 2.0 * std::exp(-2.0 / M_PI * x) * std::sin(4.0 * x);
                case 1:
                    a = 2.0 / M_PI * x - 1.0 / 2.0; 
                    return -a * a * (-a * a + 5.0 * a - 6.0) / (a * a + 1.0) - 1.0 / 2.0;
                case 2: return std::pow(std::sin(x), 3) + std::pow(std::cos(2.0 * x), 3) - 3.0 / 10.0;
                case 3:
                    b = 4.0 / M_PI * (x - 3.0 / 10.0) - 4.0;
                    return -1.0 / 500.0 * b * b * b * b * b * b + 3.0 / 100.0 * b * b * b * b -
                           27.0 / 500.0 * b * b + 3.0 / 2.0;
                default: return std::numeric_limits<double>::quiet_NaN();
            }
        },
        blockNames[1], 10, 3,
        opt::OneDimensionalSearchArea(0.0, 2.0 * M_PI), std::vector<double>{ 4.2250023 }, 1.474, 12.442132,
        std::vector<double>{ 3.170468, 4.329008, 7.999984 })
};

#endif // ONE_DIM_CONSTRAINED_TEST_PROBLEMS_H_
