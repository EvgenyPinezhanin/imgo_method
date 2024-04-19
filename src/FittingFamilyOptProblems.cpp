#include <test_opt_problems/FittingFamilyOptProblems.h>

const size_t familySizeFitting = 100;
const size_t dimensionFitting = 4;

const std::vector<double> lowerBound{ 0.01, 0.01, 0.01, 0.01 };
const std::vector<double> upperBound{ 2.0, 2.0, 2.0, 2.0 };
const opt::MultiDimensionalSearchArea searchAreaFitting(dimensionFitting, lowerBound, upperBound);

const double alphaFitting = 0.03, deltaFitting = 0.3,
             leftBoundWindowFitting = 1.0, rightBoundWindowFitting = 10.0;

// std::mt19937_64 gen(30032001);
// std::vector<double> firstValuesFitting(100), secondValuesFitting(100);
// 
// for (size_t i = 0; i < 100; ++i) {
//     firstValuesFitting[i] = ((double)(gen() - gen.min()) / (gen.max() - gen.min()) - 0.5) * 20.0;
//     secondValuesFitting[i] = ((double)(gen() - gen.min()) / (gen.max() - gen.min()) - 0.5) * 20.0;
// }

const double firstPointFitting = 13.0, secondPointFitting = 16.65,
             lastPointFitting = 18.0;
const std::vector<double> firstValuesFitting{
    -7.78107, -4.90685, -3.37658, -2.73023, -5.94084,
     7.96594, -1.50398, -7.49603, -7.45911,  2.17224,
    -1.14612, -2.40455,  1.03556,  8.00782,  4.49972,
     5.09156,  6.94298, -5.16885,  1.49283, -2.41119,
     4.15064, -1.94869, -1.83866, -1.58951,  1.09116,
    -9.23451, -9.27416, -7.83039,  5.11324,  6.33287,
     5.89405,  2.29501, -4.59979, -3.95681, -2.83353,
    -8.73609, -6.50321, -3.74151,  2.46662, -2.94344,
    -8.67077, -9.00650,  2.89405,  8.98783,  6.69010,
    -8.11935,  4.96405, -6.61755,  0.33684,  0.90863,
    -7.56295, -5.02180, -7.57752,  4.88265, -4.20966,
    -4.59785, -7.61615,  2.19408,  5.14240,  9.74696,
    -9.55962, -3.08478,  1.49403, -2.43129,  4.10099,
     2.30374, -7.77526, -6.27269,  9.42625, -8.02276,
    -1.11275,  3.03413,  5.60052,  8.05824,  6.90708,
    -9.64279, -7.01972, -3.67732,  4.85766,  9.47868,
     3.20472,  1.15996, -8.74756,  0.87256, -4.25306,
    -5.93276,  5.78546, -6.71757, -9.52325, -4.38575,
    -0.78112,  9.99732, -5.61775,  8.03561, -7.03563,
     1.69774,  3.79953, -9.02108,  5.05960, -9.34146 
};
const std::vector<double> secondValuesFitting{
    -0.00290, -4.34861,  9.67357,  2.97910, -8.32345,
     9.49470, -5.08374, -5.55973, -3.39451, -0.08234,
     9.14136, -2.34577,  4.85999,  7.25819, -2.04228,
    -9.12765, -2.78671,  0.89210, -8.20987, -4.98163,
     2.53174,  7.43376,  6.14034, -0.80006, -7.53466,
     0.59987,  0.23660,  4.67620,  8.90145, -3.70739,
    -8.46062,  5.17593, -1.67315, -0.72968, -5.54228,
    -7.86698,  6.19509, -3.52645, -1.56853, -0.59078,
     4.08909, -4.31443, -0.23901,  7.99223, -1.62053,
     5.49558,  5.35255, -6.19067, -7.54457,  0.99195,
    -9.98292, -2.90758, -3.91495,  1.03898, -5.36399,
     2.12583,  0.63012, -5.81300, -5.08053, -0.60890,
    -2.47634, -0.35167,  3.03741, -5.00285, -9.86953,
     1.45434,  0.29679, -4.39870,  1.66226,  7.20481,
    -3.55114, -3.02148, -6.26793,  0.21012, -0.00515,
    -4.42785,  3.84965, -8.90070,  4.50456, -4.65675,
     8.16895,  6.88843,  6.25763, -5.37536, -1.90303,
     9.42028, -3.51618,  1.82551, -0.52529, -9.74838,
     6.00976,  2.57032, -5.58742, -5.03271,  7.73093,
     5.69401, -1.28251, -9.38067, -0.45331, -6.41473
};

const std::vector<double> windowPointsFitting{ 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
                                               5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5 };

const std::vector<std::vector<std::vector<double>>> optimalPointsFitting {
    std::vector<std::vector<double>>{ std::vector<double>{ 0.0238464, 0.633575, 1.06063, 1.17043 } }, // 0
    std::vector<std::vector<double>>{ std::vector<double>{ 0.7861, 1.3433, 1.6617, 1.8209 } }, // 1
    std::vector<std::vector<double>>{ std::vector<double>{ 0.0299, 0.0498, 0.0896, 0.1294 } }, // 2
    std::vector<std::vector<double>>{ std::vector<double>{ 0.01, 0.0299, 0.0498, 0.1095 } }, // 3
    std::vector<std::vector<double>>{ std::vector<double>{ 0.9055, 1.5025, 1.602, 1.7413 } }, // 4
    std::vector<std::vector<double>>{ std::vector<double>{ 0.01, 0.0498, 0.0896, 1.1642 } }, // 5
    std::vector<std::vector<double>>{ std::vector<double>{ 1.6219, 1.7413, 1.8607, 1.9801 } }, // 6
    std::vector<std::vector<double>>{ std::vector<double>{ 0.2886, 0.4876, 1.3035, 1.403 } }, // 7
    std::vector<std::vector<double>>{ std::vector<double>{ 0.3085, 0.5075, 1.2836, 1.3831 } }, // 8
    std::vector<std::vector<double>>{ std::vector<double>{ 0.01, 0.0299, 0.0498, 0.1493 } }, // 9

    std::vector<std::vector<double>>{ std::vector<double>{ 0.01, 0.0299, 0.0498, 0.1692 } }, // 10
    std::vector<std::vector<double>>{ std::vector<double>{ 0.9851, 1.602, 1.7214, 1.9005 } }, // 11
    std::vector<std::vector<double>>{ std::vector<double>{ 1.6219, 1.7413, 1.8607, 1.9801 } }, // 12
    std::vector<std::vector<double>>{ std::vector<double>{ 0.3483, 0.3881, 1.2836, 1.3632 } }, // 13
    std::vector<std::vector<double>>{ std::vector<double>{ 0.01, 0.0299, 0.0697, 0.1294 } }, // 14
    std::vector<std::vector<double>>{ std::vector<double>{ 0.01, 0.8259, 0.9055, 0.9851 } }, // 15
    std::vector<std::vector<double>>{ std::vector<double>{ 0.01, 0.7463, 1.0647, 1.1642 } }, // 16
    std::vector<std::vector<double>>{ std::vector<double>{ 0.01, 0.0299, 0.0697, 0.1294 } }, // 17
    std::vector<std::vector<double>>{ std::vector<double>{ 0.01, 0.0299, 0.0498, 0.1692 } }, // 18
    std::vector<std::vector<double>>{ std::vector<double>{ 1.94024, 1.82558, 1.71093, 1.21828 } }, // 19

    std::vector<std::vector<double>>{ std::vector<double>{ 1.25132, 0.747019, 1.83044, 1.69246 } }, // 2
    std::vector<std::vector<double>>{ std::vector<double>{ 0.133889, 0.0104858, 0.196077, 0.0726733 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 0.255349, 0.0124292, 0.0134009, 0.196077 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.65165, 0.66637, 1.50007, 0.720784 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 0.0668433, 0.130974, 0.0104858, 0.197048 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 0.0182593, 0.641106, 1.1318, 1.05601 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.13861, 0.025061, 0.634304, 1.0667 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.15318, 0.998684, 0.0202026, 0.824753 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.63319, 1.75854, 1.00451, 1.50882 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.06767, 0.0104858, 1.19205, 0.756736 } },

    std::vector<std::vector<double>>{ std::vector<double>{ 1.02395, 0.935525, 1.125, 0.0308911 } }, // 3
    std::vector<std::vector<double>>{ std::vector<double>{ 1.97036, 1.83822, 1.29019, 1.72647 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 0.833499, 1.6672, 1.86348, 1.23966 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.70024, 0.815037, 1.19302, 1.8217 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.68372, 1.95579, 1.81975, 1.24355 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 0.384583, 1.29505, 0.353489, 1.43302 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.00451, 0.0104858, 0.943298, 1.0667 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.65262, 1.39221, 1.79157, 0.806292 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 0.0648999, 0.130974, 0.0124292, 0.197048 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.64582, 1.76631, 1.11529, 0.712039 } },

    std::vector<std::vector<double>>{ std::vector<double>{ 1.00646, 0.759651, 0.0143726, 1.20274 } }, // 4
    std::vector<std::vector<double>>{ std::vector<double>{ 0.508958, 1.26201, 0.285471, 1.3757 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 0.0648999, 0.130974, 0.0124292, 0.199963 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 0.0134009, 0.450657, 0.145549, 1.20177 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.22606, 0.69552, 1.03172, 0.0114575 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 0.885969, 1.06184, 0.0143726, 0.948157 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.81295, 1.43983, 0.881111, 1.62736 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 0.0785034, 0.119314, 0.159153, 1.08031 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 0.0134009, 0.0143726, 0.385554, 0.130002 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.65165, 0.993826, 1.83433, 1.94899 } },

    std::vector<std::vector<double>>{ std::vector<double>{ 1.27367, 0.0172876, 0.136804, 0.0833618 } }, // 5
    std::vector<std::vector<double>>{ std::vector<double>{ 1.70412, 0.814065, 1.27853, 1.81101 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.28047, 0.358347, 0.420535, 1.38444 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 0.887913, 1.70412, 1.84502, 1.18913 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.63222, 1.49618, 1.76631, 0.806292 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 0.0726733, 0.196077, 0.133889, 0.0104858 },
                                      std::vector<double>{ 0.0726733, 0.134861, 0.0104858, 0.197048 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.06476, 0.0153442, 0.641106, 1.13472 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 0.133889, 0.0104858, 0.196077, 0.0726733 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.12986, 0.0104858, 1.00549, 0.943298 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.14055, 0.642078, 0.0240894, 1.07059 } },

    std::vector<std::vector<double>>{ std::vector<double>{ 1.33489, 0.276726, 1.22897, 0.545881 } }, // 6
    std::vector<std::vector<double>>{ std::vector<double>{ 1.64971, 1.76145, 1.10265, 0.705237 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.98591, 1.73716, 1.86251, 1.20274 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.94024, 1.82558, 1.71093, 1.21828 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.21148, 0.952043, 1.05504, 0.625559 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.87222, 1.75757, 0.947185, 1.44371 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 0.0211743, 0.643049, 1.14249, 1.0667 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 0.478835, 1.40679, 1.31156, 0.311707 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.31059, 0.604182, 1.21245, 0.263123 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.12986, 1.00549, 0.0104858, 0.881111 } },

    std::vector<std::vector<double>>{ std::vector<double>{ 1.86542, 1.63125, 1.75659, 1.99563 } }, // 7
    std::vector<std::vector<double>>{ std::vector<double>{ 0.0726733, 0.196077, 0.133889, 0.0104858 },
                                      std::vector<double>{ 0.133889, 0.0104858, 0.196077, 0.0726733 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.12986, 0.0104858, 1.00549, 0.943298 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.14638, 0.638191, 0.0172876, 1.07253 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.14929, 0.0104858, 0.626531, 1.06865 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.26104, 0.504099, 0.282556, 1.34072 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 0.0104858, 1.12986, 1.00451, 0.817952 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.3961, 1.67594, 1.93247, 1.78672 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.34169, 0.404988, 0.370007, 1.45926 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.00937, 0.758679, 0.0240894, 1.19108 } },

    std::vector<std::vector<double>>{ std::vector<double>{ 1.41942, 1.82461, 1.93927, 1.70412 } }, // 8
    std::vector<std::vector<double>>{ std::vector<double>{ 1.73716, 1.9898, 1.60793, 1.84988 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 0.885969, 1.06087, 0.970505, 0.0279761 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 0.0134009, 0.0143726, 0.385554, 0.130002 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.24355, 1.85182, 1.67206, 0.801433 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.01229, 0.94427, 0.107654, 1.07739 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.00451, 0.881111, 0.0104858, 1.06767 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.2348, 0.70718, 1.02978, 0.0231177 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 0.025061, 0.640134, 1.13666, 1.06379 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.42914, 1.58655, 1.89652, 1.77603 },
                                      std::vector<double>{ 1.777, 1.88971, 1.44274, 1.57878} },

    std::vector<std::vector<double>>{ std::vector<double>{ 0.130974, 0.0104858, 0.0143726, 0.254377 } }, // 9
    std::vector<std::vector<double>>{ std::vector<double>{ 0.502156, 1.23772, 1.33294, 0.312678 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 0.0765601, 0.0172876, 0.0240894, 1.125 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.06767, 0.0104858, 1.19205, 0.756736 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.00451, 0.0104858, 0.943298, 1.0667 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.74008, 1.6021, 1.97717, 1.85279 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 0.127087, 0.0124292, 0.0755884, 0.194133 },
                                      std::vector<double>{ 0.0726733, 0.134861, 0.0104858, 0.197048 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 0.197048, 0.259236, 0.0104858, 1.12986 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 0.128059, 0.0104858, 0.073645, 0.194133 },
                                      std::vector<double>{ 0.0726733, 0.134861, 0.0104858, 0.197048 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 0.31365, 1.37958, 1.28047, 0.44677 } }
};

const std::vector<double> optimalValuesFitting {
     -5.15377, // 0
    -11.163,   // 1
    -21.2391,  // 2
     -9.76801, // 3
    -20.2529,  // 4
    -11.8459,  // 5
    -18.1032,  // 6
     -7.69223, // 7
     -6.46851, // 8
     -4.13564, // 9

    -19.046,   // 10
     -7.59954, // 11
    -17.9313,  // 12
     -8.90044, // 13
     -9.06334, // 14
    -11.5544,  // 15
     -7.24347, // 16
     -7.70992, // 17
    -17.7382,  // 18
    -15.6066,  // 19

     -7.30814, -15.3084, - 13.0607,   -3.12091, -14.5826, // 2
     -5.95241,  -6.00309,  -9.77678, -20.4846,   -7.81804,

    -13.5663,  -16.0088,   -5.49887,  -5.06972, -16.3472, // 3
     -9.28037, -11.4542,   -9.25541,  -5.42001,  -4.18594,

     -9.0378,   -7.3331,   -3.50313,  -8.93695,  -5.54441, // 4
    -10.5604,  -13.2762,   -7.41819, -13.8957,   -3.77639,

    -12.3072,   -8.63458,  -6.63607,  -4.7894,  -14.3287, // 5
     -8.67452,  -5.16556, -12.6835,   -9.58347,  -6.61801,

     -5.8946,   -4.46882, -10.4205,  -15.6701,  -12.4687, // 6
     -4.65406,  -5.11613,  -6.59087,  -5.65083, -13.1604,

    -12.811,    -8.60641, -11.1435,   -5.13746,  -4.60812, // 7
     -7.57139,  -8.33871, -25.7082,   -5.87832,  -9.98561,

    -24.6592,  -25.5043,  -11.9244,  -10.6042,   -5.93141, // 8
    -13.9058,   -7.58089,  -5.77671,  -5.50259, -26.3171,

    -11.5467,   -6.26045,  -7.03579, -10.1698,  -13.3914, // 9
    -19.917,    -6.33359, -10.9576,   -6.19937,  -8.9841
};

FittingFamilyOptProblems::FittingFamilyOptProblems(
    const GsaMethod<OneDimensionalSupportiveOptProblem> &_gsa, bool _isSortX)
    : BaseFittingFamilyOptProblems(familySizeFitting, dimensionFitting, searchAreaFitting,
      alphaFitting, deltaFitting, leftBoundWindowFitting, rightBoundWindowFitting,
      firstPointFitting, firstValuesFitting, secondPointFitting, secondValuesFitting,
      lastPointFitting, windowPointsFitting, _gsa, _isSortX, optimalPointsFitting, optimalValuesFitting) {};
