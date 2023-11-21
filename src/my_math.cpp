#include <my_math.h>

void mnk::solve(std::vector<double> &x) {
    std::vector<std::vector<double>> A_tmp(A);
    std::vector<double> b_tmp(b);
    size_t size =  A_tmp.size(), length;
    double coeff;

    for (size_t i = 0; i < size; ++i) { 
        for (size_t j = i; j < size; ++j) {
            if (A_tmp[j][i] != 0) {
                length = A_tmp[j].size();
                coeff = A_tmp[j][i];
                for (size_t k = i; k < length; ++k) {
                    A_tmp[j][k] /= coeff;
                }
                b_tmp[j] /= coeff;
                for (size_t k = j + 1; k < size; ++k) {
                    coeff = A_tmp[k][i];
                    for (size_t l = i; l < size; ++l) {
                        A_tmp[k][l] -= coeff * A_tmp[j][l];
                    }
                    b_tmp[k] -= coeff * b_tmp[j];
                }
                swap(A_tmp[j], A_tmp[i]);
                break;
            }
        }
    }

    for (size_t i = size - 1; i > 0; --i) {
        for (size_t j = 0; j < i; ++j) {
            b_tmp[j] -= A_tmp[j][i] * b_tmp[i];
            A_tmp[j][i] = 0.0;
        }
    }
    x = b_tmp;
}

double euclideanDistance(const std::vector<double> &val1, const std::vector<double> &val2) {
    double res = 0.0;
    size_t size = val1.size();
    for (int i = 0; i < size; i++) {
        res += (val1[i] - val2[i]) * (val1[i] - val2[i]);
    }
    return sqrt(res);
}
