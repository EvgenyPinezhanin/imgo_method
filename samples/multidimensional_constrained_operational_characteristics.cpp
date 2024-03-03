#if defined( _MSC_VER )
    #define _CRT_SECURE_NO_WARNINGS
    #define PROC_BIND
#else
    #define PROC_BIND proc_bind(master)
#endif

#include <iostream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <functional>
#include <random>

#include <opt_methods/MggsaMethod.h>
#include <opt_methods/GsaMethod.h>
#include <opt_problems/ConstrainedOptProblem.h>
#include <opt_problems/OptProblem.h>
#include <gnuplot/OutputFile.h>
#include <gnuplot/Script.h>
#include <MyMath.h>
#include <omp.h>

using namespace std;

#define CALC
#define DRAW

const int familyNumber = 0; // 0 - Fitting
const int displayType = 0; // 0 - application, 1 - png, 2 - png(notitle)

size_t N = 3;
size_t number_coefficients = 2 * N + 2;
std::vector<double> coefficients(number_coefficients);
std::vector<double> omega(N + 1);
double alpha = 0.03, delta = 0.3, a = 1.0, b = 10.0;
size_t number_window_points = 10;
std::vector<point> q{ point( a,      0.0 ),
                            point( 2.5,    0.0 ),
                            point( 4.0,    0.0 ),
                            point( 5.5,    0.0 ),
                            point( 7.0,    0.0 ),
                            point( 8.5,    0.0 ),
                            point( b,      0.0 ),
                            point( 13.0,   7.65 ),
                            point( 16.65, -9.86 ),
                            point( 18.0,   0.0  ) };

bool is_min;
const OneDimensionalOptProblem u(
    [] (double t) -> double {
        double result = 0.0;

        for (size_t i = 0; i < N + 1; ++i) {
            result += coefficients[2 * i] * std::sin(omega[i] * t) +
                      coefficients[2 * i + 1] * std::cos(omega[i] * t) ;
        }

        return is_min ? result : -result;
    },
    "Sample test problem family", 1, opt::OneDimensionalSearchArea(a, b));
double accuracy_in = 0.02, reliability_in = 2.5;
int maxTrials_in = 2000, maxFevals_in = 2000;
GsaMethod<OneDimensionalOptProblem>::Parameters parameters(accuracy_in, 0.0, maxTrials_in, maxFevals_in, reliability_in);
GsaMethod<OneDimensionalOptProblem>::Result result;

std::vector<double> optValue { 1e+13, // 0
                               -41.6923,
                               1e+13,
                               1e+13,
                               -21.0637,
                               1e+13,
                               -21.0637,
                               1e+13,
                               -97.6296,
                               -21.3581,
                               -41.124,
                               1e+13,

                               1e+13, // 1
                               -5.74494,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               -83.3182,

                               1e+13, // 2
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,

                               1e+13, // 3
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,

// 30. minVal = 1e+13 minX = (0, 0, 0, 0)
// 31. minVal = 1e+13 minX = (0, 0, 0, 0)
// 32. minVal = -85.77 minX = (0.41, 0.61, 0.51, 0.56)
// 33. minVal = 1e+13 minX = (0, 0, 0, 0)
// 34. minVal = -49.133 minX = (0.61, 0.66, 0.36, 0.41)
// 35. minVal = -31.9657 minX = (0.56, 0.51, 0.61, 0.41)
// 36. minVal = 1e+13 minX = (0, 0, 0, 0)
// 37. minVal = 1e+13 minX = (0, 0, 0, 0)
// 38. minVal = 1e+13 minX = (0, 0, 0, 0)
// 39. minVal = 1e+13 minX = (0, 0, 0, 0)
                            
                               1e+13, // 4
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,

                               1e+13, // 5
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,

                               1e+13, // 6
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,

                               1e+13, // 7
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,

                               1e+13, // 8
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,

                               1e+13, // 9
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
                               1e+13,
};

std::vector<std::vector<double>> Xopt { std::vector<double>{0.259236, 0.258264, 0.632361, 0.133889}, // 0
                                        std::vector<double>{1.00549, 0.755764, 0.382639, 0.881111},
                                        std::vector<double>{1.99951, 1.75174, 1.37764, 1.12986},
                                        std::vector<double>{0.507014, 1.72744, 1.45246, 0.305876},
                                        std::vector<double>{0.01, 1.01, 0.76, 0.66},
                                        std::vector<double>{1.00549, 1.00549, 1.37764, 0.133889},
                                        std::vector<double>{0.61, 0.41, 0.46, 0.56},
                                        std::vector<double>{1.75076, 0.259236, 0.631389, 1.37861},
                                        std::vector<double>{0.507986, 0.507986, 0.881111, 1.37861},
                                        std::vector<double>{1.25424, 1.50201, 1.00549, 0.755764},

//1 std::vector<double>{0.259236, 0.258264, 0.632361, 0.133889}
//2 std::vector<double>{1.00549, 0.755764, 0.382639, 0.881111}
//2 std::vector<double>{1.99951, 1.75174, 1.37764, 1.12986}
//1 std::vector<double>{0.507986, 1.75174, 1.50299, 0.259236}
//1 std::vector<double>{0.0104858, 1.00549, 0.756736, 0.756736}
//1 std::vector<double>{1.00451, 1.00549, 1.37861, 0.133889}
//1 std::vector<double>{0.563372, 0.404988, 0.524504, 0.522561}
//2 std::vector<double>{1.75076, 0.259236, 0.631389, 1.37861}
//2 std::vector<double>{0.507986, 0.507986, 0.881111, 1.37861}
//2 std::vector<double>{1.25424, 1.50201, 1.00549, 0.755764}

                                        std::vector<double>{0.285471, 0.0182593, 1.47384, 1.4991}, // 1
                                        std::vector<double>{0.51, 0.71, 1.06, 0.06},
                                        std::vector<double>{0.0308911, 0.515759, 0.509929, 0.756736},
                                        std::vector<double>{1.26687, 0.115427, 1.14444, 1.17747},
                                        std::vector<double>{0.630417, 0.879167, 1.25326, 0.755764},
                                        std::vector<double>{1.62542, 1.74979, 0.754792, 1.12986},
                                        std::vector<double>{1.12986, 0.755764, 0.756736, 0.880139},
                                        std::vector<double>{0.877224, 1.00451, 1.27173, 0.636248},
                                        std::vector<double>{1.75174, 0.383611, 0.507986, 1.62736},
                                        std::vector<double>{0.46, 0.61, 0.41, 0.56},

// 2 std::vector<double>{0.285471, 0.0182593, 1.47384, 1.4991}
// 2 std::vector<double>{1.43205, 1.76534, 0.504099, 1.35335}
// 2 std::vector<double>{0.0308911, 0.515759, 0.509929, 0.756736}
// 2 std::vector<double>{1.26687, 0.115427, 1.14444, 1.17747}
// 2 std::vector<double>{0.630417, 0.879167, 1.25326, 0.755764}
// 2 std::vector<double>{1.62542, 1.74979, 0.754792, 1.12986}
// 2 std::vector<double>{1.12986, 0.755764, 0.756736, 0.880139}
// 2 std::vector<double>{0.877224, 1.00451, 1.27173, 0.636248}
// 2 std::vector<double>{1.75174, 0.383611, 0.507986, 1.62736}
// 2 std::vector<double>{1.12889, 1.37764, 0.258264, 1.50201}

                                        std::vector<double>{1.4991, 1.99757, 0.754792, 0.756736}, // 2
                                        std::vector<double>{1.00549, 0.507986, 0.629446, 0.132917},
                                        std::vector<double>{1.25326, 1.87611, 1.00549, 1.87514},
                                        std::vector<double>{0.0308911, 0.987996, 1.98105, 1.75854},
                                        std::vector<double>{1.25326, 1.37861, 1.50299, 0.382639},
                                        std::vector<double>{0.133889, 0.0104858, 0.258264, 0.0104858},
                                        std::vector<double>{1.28533, 1.80129, 1.7051, 1.41748},
                                        std::vector<double>{0.507986, 1.50201, 1.50299, 0.258264},
                                        std::vector<double>{1.00549, 1.25424, 1.87514, 0.382639},
                                        std::vector<double>{1.25326, 1.37764, 0.258264, 0.880139},

// 1 std::vector<double>{1.74979, 1.50493, 1.00257, 1.74979}
// 2 std::vector<double>{1.00549, 0.507986, 0.629446, 0.132917}
// 2 std::vector<double>{1.25326, 1.87611, 1.00549, 1.87514}
// 2 std::vector<double>{0.0308911, 0.987996, 1.98105, 1.75854}
// 1 std::vector<double>{1.72842, 1.52436, 0.246604, 1.05407}
// 1 std::vector<double>{0.133889, 0.0104858, 0.258264, 0.0104858}
// 2 std::vector<double>{0.500212, 0.501184, 1.02006, 1.75465}
// 1 std::vector<double>{0.507986, 1.50396, 1.50299, 0.257292}
// 2 std::vector<double>{1.00549, 1.25424, 1.87514, 0.382639}
// 1 std::vector<double>{1.25326, 1.37764, 0.258264, 0.880139}

                                        std::vector<double>{1.42622, 1.00063, 1.05213, 1.04144}, // 3
                                        std::vector<double>{1.50201, 1.99951, 1.00451, 0.755764},
                                        std::vector<double>{1.00549, 0.259236, 0.881111, 1.12986},
                                        std::vector<double>{1.50201, 0.258264, 0.881111, 1.37861},
                                        std::vector<double>{0.383611, 1.50299, 0.756736, 1.12986},
                                        std::vector<double>{1.00549, 0.756736, 0.880139, 1.87514},
                                        std::vector<double>{0.262151, 0.152351, 0.488552, 0.356404},
                                        std::vector<double>{0.717869, 1.45052, 1.52242, 1.29893},
                                        std::vector<double>{0.507986, 0.880139, 0.755764, 0.632361},
                                        std::vector<double>{1.5127, 1.49133, 1.48938, 0.265066},

// 1 std::vector<double>{1.50201, 1.00451, 1.00451, 1.00451}
// 1 std::vector<double>{1.50201, 1.99951, 1.00451, 0.755764}
// 2 std::vector<double>{1.00549, 0.259236, 0.881111, 1.12986}
// 1 std::vector<double>{1.50201, 0.258264, 0.881111, 1.37861}
// 2 std::vector<double>{0.383611, 1.50299, 0.756736, 1.12986}
// 2 std::vector<double>{1.00549, 0.756736, 0.880139, 1.87514}
// 1 std::vector<double>{0.259236, 0.133889, 0.507986, 0.382639}
// 1 std::vector<double>{0.726614, 1.4408, 1.53991, 1.26007}
// 1 std::vector<double>{0.507986, 0.880139, 0.755764, 0.632361}
// 2 std::vector<double>{1.5127, 1.49133, 1.48938, 0.265066}

                                        std::vector<double>{0.0104858, 1.00451, 1.25326, 1.75174}, // 4
                                        std::vector<double>{0.507014, 1.00549, 1.00451, 1.75174},
                                        std::vector<double>{0.520618, 0.522561, 1.14055, 1.61376},
                                        std::vector<double>{0.258264, 0.631389, 0.507014, 1.37764},
                                        std::vector<double>{0.0211743, 1.75854, 1.37084, 1.61959},
                                        std::vector<double>{0.880139, 1.25326, 1.50299, 0.632361},
                                        std::vector<double>{1.99951, 1.99951, 0.630417, 1.13083},
                                        std::vector<double>{0.51, 0.71, 0.41, 0.36},
                                        std::vector<double>{0.134861, 0.507986, 1.25424, 0.632361},
                                        std::vector<double>{0.564343, 1.53894, 1.96356, 1.7187},

// 2 std::vector<double>{0.0104858, 1.00451, 1.25326, 1.75174}
// 2 std::vector<double>{0.507014, 1.00549, 1.00451, 1.75174}
// 1 std::vector<double>{0.510901, 0.520618, 1.1318, 1.62153}
// 1 std::vector<double>{0.254377, 0.633333, 0.511873, 1.3757}
// 2 std::vector<double>{0.900544, 1.63708, 1.00451, 0.128059}
// 2 std::vector<double>{1.78672, 0.335027, 1.21926, 1.68275}
// 1 std::vector<double>{1.99951, 1.99854, 0.631389, 1.13083}
// 1 std::vector<double>{0.507986, 0.756736, 0.383611, 0.383611}
// 2 std::vector<double>{0.134861, 0.507986, 1.25424, 0.632361}
// 2 std::vector<double>{0.564343, 1.53894, 1.96356, 1.7187}


//                                         40. minVal = 1e+13 minX = (0, 0, 0, 0)
// 41. minVal = -31.5178 minX = (0.81, 0.96, 0.76, 0.71)
// 42. minVal = 1e+13 minX = (0, 0, 0, 0)
// 43. minVal = 1e+13 minX = (0, 0, 0, 0)
// 44. minVal = 1e+13 minX = (0, 0, 0, 0)
// 45. minVal = 1e+13 minX = (0, 0, 0, 0)
// 46. minVal = 1e+13 minX = (0, 0, 0, 0)
// 47. minVal = -55.4893 minX = (0.51, 0.71, 0.41, 0.36)
// 48. minVal = 1e+13 minX = (0, 0, 0, 0)
// 49. minVal = 1e+13 minX = (0, 0, 0, 0)

                                        std::vector<double>{0.755764, 1.50299, 1.12889, 1.37861}, // 5
                                        std::vector<double>{1.75174, 0.507986, 1.99951, 1.75076}, 
                                        std::vector<double>{1.90137, 0.522561, 1.80226, 0.490496}, 
                                        std::vector<double>{1.83239, 0.780056, 1.48064, 0.744104}, 
                                        std::vector<double>{0.631389, 0.507986, 0.756736, 1.37764}, 
                                        std::vector<double>{0.259236, 0.880139, 1.00549, 0.880139},
                                        std::vector<double>{0.371951, 0.745076, 0.0114575, 0.375837},
                                        std::vector<double>{1.00451, 1.00549, 1.00549, 0.258264},
                                        std::vector<double>{1.00451, 0.756736, 1.00451, 0.756736},
                                        std::vector<double>{0.134861, 1.50201, 1.00451, 0.383611},

// 2 std::vector<double>{0.755764, 1.50299, 1.12889, 1.37861}
// 2 std::vector<double>{1.75174, 0.507986, 1.99951, 1.75076}
// 2 std::vector<double>{1.90137, 0.522561, 1.80226, 0.490496}
// 2 std::vector<double>{1.83239, 0.780056, 1.48064, 0.744104}
// 2 std::vector<double>{0.631389, 0.507986, 0.756736, 1.37764}
// 2 std::vector<double>{0.259236, 0.880139, 1.00549, 0.880139}
// 2 std::vector<double>{0.371951, 0.745076, 0.0114575, 0.375837}
// 2 std::vector<double>{1.00451, 1.00549, 1.00549, 0.258264}
// 2 std::vector<double>{1.00451, 0.756736, 1.00451, 0.756736}
// 2 std::vector<double>{0.134861, 1.50201, 1.00451, 0.383611}

//                                         50. minVal = -31.1678 minX = (0.96, 0.71, 0.66, 0.01)
// 51. minVal = -82.29 minX = (0.46, 0.41, 0.51, 0.66)
// 52. minVal = -36.5842 minX = (1.36, 1.26, 1.06, 1.16)
// 53. minVal = 1e+13 minX = (0, 0, 0, 0)
// 54. minVal = -69.3791 minX = (0.36, 0.41, 0.66, 0.61)
// 55. minVal = 1e+13 minX = (0, 0, 0, 0)
// 56. minVal = 1e+13 minX = (0, 0, 0, 0)
// 57. minVal = 1e+13 minX = (0, 0, 0, 0)
// 58. minVal = 1e+13 minX = (0, 0, 0, 0)
// 59. minVal = 1e+13 minX = (0, 0, 0, 0)

                                        std::vector<double>{0.66, 0.46, 0.31, 0.36}, // 6
                                        std::vector<double>{0.756736, 1.99951, 1.50201, 1.25424},
                                        std::vector<double>{1.25326, 1.12986, 1.25424, 1.37764},
                                        std::vector<double>{0.56, 0.41, 0.46, 0.61},
                                        std::vector<double>{1.41, 0.51, 0.06, 0.96},
                                        std::vector<double>{1.37764, 1.00549, 0.756736, 1.37764},
                                        std::vector<double>{1.12986, 1.75174, 1.00549, 1.87611},
                                        std::vector<double>{0.66, 0.56, 0.36, 0.41},
                                        std::vector<double>{0.758679, 0.630417, 1.25326, 1.87028},
                                        std::vector<double>{1.14152, 1.2999, 1.86445, 1.28727},

// 2 std::vector<double>{0.292273, 0.477864, 0.558513, 0.727585}
// 2 std::vector<double>{0.756736, 1.99951, 1.50201, 1.25424}
// 2 std::vector<double>{1.25326, 1.12986, 1.25424, 1.37764}
// 2 std::vector<double>{0.120286, 0.264094, 0.764509, 0.137776}
// 2 std::vector<double>{1.00549, 1.25424, 0.880139, 0.880139}
// 2 std::vector<double>{1.37764, 1.00549, 0.756736, 1.37764}
// 2 std::vector<double>{1.12986, 1.75174, 1.00549, 1.87611}
// 2 std::vector<double>{0.821838, 0.0435229, 1.08808, 0.946213}
// 2 std::vector<double>{0.758679, 0.630417, 1.25326, 1.87028}
// 2 std::vector<double>{1.14152, 1.2999, 1.86445, 1.28727}

//                                         60. minVal = -36.0963 minX = (0.66, 0.46, 0.31, 0.36)
// 61. minVal = 1e+13 minX = (0, 0, 0, 0)
// 62. minVal = 1e+13 minX = (0, 0, 0, 0)
// 63. minVal = -81.5261 minX = (0.56, 0.41, 0.46, 0.61)
// 64. minVal = -0.0476254 minX = (1.41, 0.51, 0.06, 0.96)
// 65. minVal = 1e+13 minX = (0, 0, 0, 0)
// 66. minVal = 1e+13 minX = (0, 0, 0, 0)
// 67. minVal = -36.476 minX = (0.66, 0.56, 0.36, 0.41)
// 68. minVal = 1e+13 minX = (0, 0, 0, 0)
// 69. minVal = 1e+13 minX = (0, 0, 0, 0)

                                        std::vector<double>{0.91, 1.01, 0.76, 0.81}, // 7
                                        std::vector<double>{0.504099, 0.635276, 1.25132, 1.63028},
                                        std::vector<double>{1.75174, 0.134861, 0.0104858, 0.632361},
                                        std::vector<double>{0.0104858, 0.507986, 1.50201, 1.75076},
                                        std::vector<double>{1.51076, 1.97911, 0.495354, 0.272839},
                                        std::vector<double>{0.259236, 1.25326, 0.632361, 1.12889},
                                        std::vector<double>{0.511873, 1.00451, 1.25132, 1.25132},
                                        std::vector<double>{1.50299, 0.259236, 0.632361, 0.632361},
                                        std::vector<double>{1.89069, 1.15027, 0.2845, 0.500212},
                                        std::vector<double>{1.48355, 0.148464, 1.26492, 1.29019},

// 1 std::vector<double>{0.991882, 1.01715, 0.746047, 0.843215}
// 2 std::vector<double>{0.504099, 0.635276, 1.25132, 1.63028}
// 2 std::vector<double>{1.75174, 0.134861, 0.0104858, 0.632361}
// 2 std::vector<double>{0.0104858, 0.507986, 1.50201, 1.75076}
// 2 std::vector<double>{0.881111, 0.507986, 1.00451, 1.12889}
// 2 std::vector<double>{1.00549, 0.383611, 1.75076, 1.37764}
// 2 std::vector<double>{1.00451, 1.50201, 1.25326, 1.25326}
// 2 std::vector<double>{0.382639, 1.25521, 0.507014, 1.87708}
// 2 std::vector<double>{1.70315, 1.62153, 1.78186, 1.75465}
// 2 std::vector<double>{0.635276, 0.756736, 1.25229, 0.380696}

// 2 std::vector<double>{1.51076, 1.97911, 0.495354, 0.272839}
// 2 std::vector<double>{0.259236, 1.25326, 0.632361, 1.12889}
// 2 std::vector<double>{0.511873, 1.00451, 1.25132, 1.25132}
// 2 std::vector<double>{1.50299, 0.259236, 0.632361, 0.632361}
// 2 std::vector<double>{1.89069, 1.15027, 0.2845, 0.500212}
// 2 std::vector<double>{1.48355, 0.148464, 1.26492, 1.29019}

// 70. minVal = -39.347 minX = (0.91, 1.01, 0.76, 0.81)
// 71. minVal = 1e+13 minX = (0, 0, 0, 0)
// 72. minVal = 1e+13 minX = (0, 0, 0, 0)
// 73. minVal = 1e+13 minX = (0, 0, 0, 0)
// 74. minVal = 1e+13 minX = (0, 0, 0, 0)
// 75. minVal = -36.9009 minX = (0.61, 0.91, 0.66, 0.51)
// 76. minVal = 1e+13 minX = (0, 0, 0, 0)
// 77. minVal = -32.0589 minX = (0.71, 0.61, 0.41, 0.56)
// 78. minVal = 1e+13 minX = (0, 0, 0, 0)
// 79. minVal = 1e+13 minX = (0, 0, 0, 0)

                                        std::vector<double>{0.501184, 0.503127, 1.51756, 1.77117}, // 8
                                        std::vector<double>{1.36404, 0.605154, 1.48647, 0.966619},
                                        std::vector<double>{1.00549, 1.50201, 0.259236, 1.25326},
                                        std::vector<double>{0.613899, 0.760623, 1.84405, 1.88291},
                                        std::vector<double>{0.76, 0.46, 0.51, 0.41},
                                        std::vector<double>{1.50201, 1.99951, 1.25424, 1.75174},
                                        std::vector<double>{0.507986, 0.507986, 1.37861, 1.37861},
                                        std::vector<double>{1.20662, 0.450657, 0.47009, 0.267981},
                                        std::vector<double>{0.56, 0.31, 0.26, 0.61},
                                        std::vector<double>{0.46, 0.51, 0.71, 0.56},

// 2 std::vector<double>{0.501184, 0.503127, 1.51756, 1.77117}
// 2 std::vector<double>{1.36404, 0.605154, 1.48647, 0.966619}
// 2 std::vector<double>{1.00549, 1.50201, 0.259236, 1.25326}
// 2 std::vector<double>{0.613899, 0.760623, 1.84405, 1.88291}
// 2 std::vector<double>{1.00549, 1.12889, 1.75076, 0.134861}
// 2 std::vector<double>{1.50201, 1.99951, 1.25424, 1.75174}
// 2 std::vector<double>{0.507986, 0.507986, 1.37861, 1.37861}
// 2 std::vector<double>{1.20662, 0.450657, 0.47009, 0.267981}
// 2 std::vector<double>{1.00549, 0.507014, 0.632361, 0.133889}
// 2 std::vector<double>{1.11529, 1.24646, 1.69538, 1.65457}

// 80. minVal = 1e+13 minX = (0, 0, 0, 0)
// 81. minVal = 1e+13 minX = (0, 0, 0, 0)
// 82. minVal = 1e+13 minX = (0, 0, 0, 0)
// 83. minVal = 1e+13 minX = (0, 0, 0, 0)
// 84. minVal = -76.3202 minX = (0.76, 0.46, 0.51, 0.41)
// 85. minVal = 1e+13 minX = (0, 0, 0, 0)
// 86. minVal = 1e+13 minX = (0, 0, 0, 0)
// 87. minVal = 1e+13 minX = (0, 0, 0, 0)
// 88. minVal = -72.0493 minX = (0.56, 0.31, 0.26, 0.61)
// 89. minVal = -28.5745 minX = (0.46, 0.51, 0.71, 0.56)

                                        std::vector<double>{1.25326, 0.507986, 0.383611, 0.133889}, // 9
                                        std::vector<double>{0.243689, 0.997712, 0.286443, 0.506042},
                                        std::vector<double>{0.01, 0.46, 0.91, 0.61},
                                        std::vector<double>{0.81, 0.51, 0.96, 0.36},
                                        std::vector<double>{1.99951, 1.87417, 0.508958, 1.12792},
                                        std::vector<double>{1.50201, 1.25326, 1.50299, 0.259236},
                                        std::vector<double>{1.01617, 1.50104, 0.482722, 0.266038},
                                        std::vector<double>{0.46, 0.66, 0.61, 0.41},
                                        std::vector<double>{1.53894, 1.1901, 1.72647, 1.28824},
                                        std::vector<double>{0.91, 0.06, 0.66, 0.76},

// 2 std::vector<double>{1.25326, 0.507986, 0.383611, 0.133889}
// 2 std::vector<double>{0.243689, 0.997712, 0.286443, 0.506042}
// 2 std::vector<double>{0.262151, 1.25715, 1.87222, 0.132917}
// 2 std::vector<double>{0.259236, 0.635276, 1.25132, 0.136804}
// 2 std::vector<double>{1.99951, 1.87417, 0.508958, 1.12792}
// 2 std::vector<double>{1.50201, 1.25326, 1.50299, 0.259236}
// 2 std::vector<double>{1.01617, 1.50104, 0.482722, 0.266038}
// 2 std::vector<double>{1.31059, 0.839329, 1.228, 1.21148}
// 2 std::vector<double>{1.53894, 1.1901, 1.72647, 1.28824}
// 2 std::vector<double>{0.639163, 1.25909, 1.01715, 0.142634}

// 90. minVal = 1e+13 minX = (0, 0, 0, 0)
// 91. minVal = 1e+13 minX = (0, 0, 0, 0)
// 92. minVal = -7.38392 minX = (0.01, 0.46, 0.91, 0.61)
// 93. minVal = -24.7469 minX = (0.81, 0.51, 0.96, 0.36)
// 94. minVal = 1e+13 minX = (0, 0, 0, 0)
// 95. minVal = 1e+13 minX = (0, 0, 0, 0)
// 96. minVal = 1e+13 minX = (0, 0, 0, 0)
// 97. minVal = -55.812 minX = (0.46, 0.66, 0.61, 0.41)
// 98. minVal = 1e+13 minX = (0, 0, 0, 0)
// 99. minVal = -34.9524 minX = (0.91, 0.06, 0.66, 0.76)
};

double u_der(std::vector<double> x) {
    double result = 0.0;
    for (size_t i = 0; i < N + 1; ++i) {
        result += coefficients[2 * i] * x[i] * std::cos(x[i] * q[number_window_points - 1].x[0]) -
                  coefficients[2 * i + 1] * x[i] * std::sin(x[i] * q[number_window_points - 1].x[0]) ;
    }
    return result;
}

std::vector<double> x_old;

// const MultiDimensionalConstrainedProblem sample_test_problem_family(
//     [] (std::vector<double> x, int index) -> double {
double problem(std::vector<double> x, int index) {
        std::vector<std::vector<double>> A(number_coefficients,
                                           std::vector<double>(number_coefficients, 0));
        std::vector<double> b(number_coefficients, 0);

        std::vector<std::function<double(double, double)>> functions{
            [] (double t, double x) -> double { return std::sin(x * t); },
            [] (double t, double x) -> double { return std::cos(x * t); }
        };

        x_old = x;

        mnk minimizer;

        GsaMethod<OneDimensionalOptProblem> gsa(u, parameters);
        double min_value, max_value;

        switch (index) {
            case 0: case 1: case 2:
                std::sort(x.begin(), x.end());
                return x[index] * (1.0 + alpha) - x[index + 1] * (1.0 - alpha);

            case 3:
                for (size_t i = 0; i < number_coefficients; ++i) {
                    for (size_t j = 0; j < number_coefficients; ++j) {
                        for (size_t k = 0; k < number_window_points; ++k) {
                            A[i][j] += functions[j % 2](q[k].x[0], x[j / 2]) * functions[i % 2](q[k].x[0], x[i / 2]);
                        }
                    }
                    for (size_t j = 0; j < number_window_points; ++j) {
                        b[i] += q[j].x[1] * functions[i % 2](q[j].x[0], x[i / 2]);
                    }
                }
                
                minimizer.set_A(A);
                minimizer.set_b(b);
                minimizer.solve(coefficients);

            case 4: case 5:

                omega = x;
                is_min = true;

                return std::abs(u.computeObjectiveFunction(q[index - 3].x[0]) - q[index - 3].x[1]) - delta;
            
            case 6:
                omega = x;

                is_min = true;
                gsa.setProblem(u);
                gsa.solve(result);
                min_value = result.value;

                is_min = false;
                gsa.setProblem(u);
                gsa.solve(result);
                max_value = -result.value;

                // std::cout << max_value - min_value - 2 * delta << "\n";

                return max_value - min_value - 2 * delta;

            case 7:
                return -std::abs(u_der(x));

            default: return std::numeric_limits<double>::quiet_NaN();
        }
    }
//    opt::MultiDimensionalSearchArea(N + 1, std::vector<double> { 0.01, 0.01, 0.01, 0.01 },
//                                           std::vector<double> { 2.0,  2.0,  2.0,  2.0 }), 2 * N + 1);
std::vector<double> A{ 0.01, 0.01, 0.01, 0.01 };
std::vector<double> B{ 2.0,  2.0,  2.0,  2.0 };

int main() {
    OutputFile vars_file("output_data/mggsa_operational_characteristics/vars.txt");
    if (!vars_file.isOpen()) std::cerr << "vars_file opening error\n";

    const int chunk = 1;

    const int numberFamily = 1;

    vector<vector<int>> K{ { 0, 5500, 25 } };

    vector<vector<double>> r{ { 3.0, 3.1, 3.2, 2.9 } };
    vector<int> key{ 1, 1, 1, 1 };
    vector<double> d{ 0.01, 0.01, 0.01, 0.01 };

    std::function<double(std::vector<double>, int)> fitting_problem_family = problem;

    // std::random_device dev;
    std::mt19937_64 gen(30032001);
    std::vector<double> point1(100), point2(100);

    for (size_t i = 0; i < 100; ++i) {
        point1[i] = ((double)(gen() - gen.min()) / (gen.max() - gen.min()) - 0.5) * 20.0;
        point2[i] = ((double)(gen() - gen.min()) / (gen.max() - gen.min()) - 0.5) * 20.0;
    }

#if defined(CALC)
    ofstream ofstr("output_data/multidimensional_constrained_operational_characteristics/operational_characteristics.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";

    int den = 11, incr = 30;
    double eps = 0.1;
    int maxFevals = 1000000;

    MggsaMethod mggsa(nullptr, -1, -1, vector<double>{}, vector<double>{}, -1.0, -1.0, den, -1, eps, -1, maxFevals, incr);

    int sizeR = r[0].size();
    vector<vector<vector<double>>> successRate(numberFamily);
    for (int i = 0; i < numberFamily; i++) {
        successRate[i].resize(sizeR);
        for (int j = 0; j < sizeR; j++) {
            successRate[i][j].resize((K[i][1] - K[i][0]) / K[i][2] + 1);
        }
    }

    // double step = 0.05;
// 
    // size_t numberFunctions = 100;
    // for (int num = 0; num < numberFunctions; num++) {
    //     q[7].x[0] = point1[num];
    //     q[8].x[0] = point2[num];
// 
    //     double minValue = 10000000000000.0;
    //     std::vector<double> Xmin(4);
    //     double tmpValue;
    //     bool b;
// 
    //     std::vector<double> X(4);
// 
    //     for (double i = A[0]; i <= B[0]; i+= step) {
    //         for (double j = A[1]; j <= B[1]; j+= step) {
    //             for (double k = A[2]; k <= B[2]; k+= step) {
    //                 for (double l = A[3]; l <= B[3]; l+= step) {
    //                     X = std::vector<double>{i, j, k, l};
    //                     b = true;
    //                     
    //                     for (int constr = 0; constr < 7; constr++) {
    //                         if (problem(X, constr) > 0.0) {
    //                             b = false;
    //                             break;
    //                         }
    //                     } 
// 
    //                     if (b) {
    //                         // std::cout << problem(X, 7) << "\n";
    //                         tmpValue = problem(X, 7);
    //                         if (tmpValue < minValue) {
    //                             Xmin = X;
    //                             minValue = tmpValue;
    //                         }
    //                     }
    //                 } 
    //             } 
    //         } 
    //     }
// 
    //     std::cout << num << ". minVal = " << minValue << " minX = (" << Xmin[0] << ", " << Xmin[1] << ", " << Xmin[2] << ", " <<  Xmin[3] << ")\n";
    // }



    double totalStartTime = omp_get_wtime();
#pragma omp parallel for schedule(dynamic, chunk) PROC_BIND num_threads(1) collapse(2) \
        shared(numberFamily, fitting_problem_family, r, key, d, sizeR, successRate, A, B) \
        firstprivate(mggsa)
    for (int i = 0; i < numberFamily; i++) {
        for (int j = 0; j < sizeR; j++) {
            
            mggsa.setNumberConstraints(2 * N + 1);
            mggsa.setMaxTrials(K[i][1]);
            mggsa.setAB(A, B);
            mggsa.setD(d[i]);
            mggsa.setKey(key[j]);
            mggsa.setR(r[i][j]);
            mggsa.setN(N + 1);

            vector<double> XOpt_o;
            int numberSuccessful, numberTrials, numberFevals;
            int numberFunctions = 100;
            vector<int> numberTrialsArray(numberFunctions);

            double startTime = omp_get_wtime();
            for (int k = 0; k < numberFunctions; k++) {
                XOpt_o = Xopt[k];
                q[7].x[0] = point1[k];
                q[8].x[0] = point2[k];
                mggsa.setF(fitting_problem_family);

                if (mggsa.solveTest(XOpt_o, numberTrials, numberFevals)) {
                    numberTrialsArray[k] = numberTrials;
                } else {
                    numberTrialsArray[k] = K[i][1] + 1;
                }
            }

            for (int k = K[i][0]; k <= K[i][1]; k += K[i][2]) {
                numberSuccessful = (int)count_if(numberTrialsArray.begin(), numberTrialsArray.end(), [k](double elem){ return elem <= k; });
                successRate[i][j][k / K[i][2]] = (double)numberSuccessful / numberFunctions;
            }
            double endTime = omp_get_wtime();
            double workTime = endTime - startTime;

            string strOutput = "Fitting model Family r = " + to_string(r[i][j]) + " key = " + to_string(key[j]) + 
                               " time: " + to_string(workTime) + " thread number: " + to_string(omp_get_thread_num()) + "\n";
            cout << strOutput;
        }
    }

    for (int i = 0; i < numberFamily; i++) {
        for (int j = 0; j < sizeR; j++) {
            for (int k = K[i][0]; k <= K[i][1]; k += K[i][2]) {
                ofstr << k << " " << successRate[i][j][k / K[i][2]] << endl;
            }
            ofstr << endl << endl;
        }
    }
    ofstr.close();
    double totalEndTime = omp_get_wtime();
    double totalWorkTime = totalEndTime - totalStartTime;
    cout << "Total time: " << totalWorkTime << endl;
#endif

    int sizeKey = key.size();
    vars_file.setVariable("numberKey", sizeKey, false);
    vars_file.initArray("familyName", numberFamily);
    vars_file.initArray("r", sizeKey * numberFamily);
    vars_file.initArray("key", sizeKey);
    for (int i = 0; i < numberFamily; i++) {
        // vars_file.setValueInArray("familyName", i + 1, problems[i].shortName);
        vars_file.setValueInArray("familyName", i + 1, "Fitting model Family");
        for (int j = 0; j < r[i].size(); j++) {
            vars_file.setValueInArray("r", (i * sizeKey) + j + 1, r[i][j]);
        }
    }
    for (int i = 0; i < sizeKey; i++) {
        vars_file.setValueInArray("key", i + 1, key[i]);
    }
    vars_file.close();

#if defined( DRAW )
    Script script("scripts/mggsa_operational_characteristics.gp");
    script.addArgs(std::vector<int>{ displayType, familyNumber });
    script.start();
    if (script.isError() == 2) std::cerr << "Error gnuplot\n";
    if (script.isError() == 1) std::cerr << "Error chmod\n";
#endif

#if defined( _MSC_VER )
    cin.get();
#endif

	return 0;
}
