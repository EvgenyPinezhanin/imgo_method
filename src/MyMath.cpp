#include <MyMath.h>

void mnk::solve(std::vector<double> &x) const {
    std::vector<std::vector<double>> ATmp(A);
    std::vector<double> BTmp(B);
    size_t size =  ATmp.size(), length;
    double coeff;

    for (size_t i = 0; i < size; ++i) { 
        for (size_t j = i; j < size; ++j) {
            if (ATmp[j][i] != 0) {
                length = ATmp[j].size();
                coeff = ATmp[j][i];
                for (size_t k = i; k < length; ++k) {
                    ATmp[j][k] /= coeff;
                }
                BTmp[j] /= coeff;
                for (size_t k = j + 1; k < size; ++k) {
                    coeff = ATmp[k][i];
                    for (size_t l = i; l < size; ++l) {
                        ATmp[k][l] -= coeff * ATmp[j][l];
                    }
                    BTmp[k] -= coeff * BTmp[j];
                }
                swap(ATmp[j], ATmp[i]);
                break;
            }
        }
    }

    for (size_t i = size - 1; i > 0; --i) {
        for (size_t j = 0; j < i; ++j) {
            BTmp[j] -= ATmp[j][i] * BTmp[i];
            ATmp[j][i] = 0.0;
        }
    }
    x = BTmp;
}

double euclideanDistance(const std::vector<double> &firstValue, const std::vector<double> &secondValue) {
    double res = 0.0;
    size_t size = firstValue.size();
    for (int i = 0; i < size; i++) {
        res += (firstValue[i] - secondValue[i]) * (firstValue[i] - secondValue[i]);
    }
    return sqrt(res);
}

double euclideanDistance(double firstValue, double secondValue) {
    return std::abs(firstValue - secondValue);
}
