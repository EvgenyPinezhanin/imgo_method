#include <gnuplot/output_file.h>

void output_file::init_array(const std::string &name, size_t number_values) {
    std::stringstream str_array;
    str_array << "array " << name << "[" << number_values << "]\n";
    stream << str_array.str();
}

void output_file::add_point(double x, double f, bool space) {
    std::stringstream str_point;
    str_point << x << " " << f << "\n";
    if (space) str_point << "\n\n";
    stream << str_point.str();
}

void output_file::add_point(const std::vector<double> &x, double f, bool space) {
    std::stringstream str_point;
    size_t dim_x = x.size();
    for (size_t i = 0; i < dim_x; ++i) {
        str_point << x[i] << " ";
    }
    str_point << f << "\n";
    if (space) str_point << "\n\n";
    stream << str_point.str();
}

void output_file::add_points(const std::vector<Trial> &trials, bool space) {
    std::stringstream str_points;
    size_t number_points = trials.size();
    for (size_t i = 0; i < number_points; ++i) {
        str_points << trials[i].x << " " << trials[i].z << "\n";
    }
    if (space) str_points << "\n\n";
    stream << str_points.str();
}

void output_file::add_points(const std::vector<opt::IndexTrial> &trials, bool space) {
    std::stringstream str_points;
    size_t number_points = trials.size();
    for (size_t i = 0; i < number_points; ++i) {
        str_points << trials[i].x << " " << trials[i].z << trials[i].nu << "\n";
    }
    if (space) str_points << "\n\n";
    stream << str_points.str();
}

void output_file::add_points(const std::vector<double> &x, double f, bool space) {
    std::stringstream str_points;
    size_t number_points = x.size();
    for (size_t i = 0; i < number_points; ++i) {
        str_points << x[i] << " " << f << "\n";
    }
    if (space) str_points << "\n\n";
    stream << str_points.str();
}

/*

void addPointGnuplot(ofstream &ofstr, const vector<double> &X) {
    size_t dimensionX = X.size();
    for (int i = 0; i < dimensionX - 1; i++) {
        ofstr << X[i] << " ";
    }
    ofstr << X[dimensionX - 1] << "\n";
    ofstr << "\n" << endl;
}

void addPointGnuplot(ofstream &ofstr, Trial trial) {
    ofstr << trial.x << " " << trial.z << "\n";
    ofstr << "\n" << endl;
}

void addPointsGnuplot(ofstream &ofstr, const vector<double> &X, const vector<double> &f) {
    size_t numberPoints = X.size();
    for (int i = 0; i < numberPoints; i++) {
        ofstr << X[i] << " " << f[i] << "\n";
    }
    ofstr << "\n" << endl;
}

void addPointsGnuplot(ofstream &ofstr, const vector<vector<double>> &X, const vector<double> &f) {
    size_t sizeX = X.size(), dimensionX = X[0].size();
    for (int i = 0; i < sizeX; i++) {
        for (int j = 0; j < dimensionX; j++) {
            ofstr << X[i][j] << " ";
        }
        ofstr << f[i] << "\n";
    }
    ofstr << "\n" << endl;
}

void addPointsGnuplot(ofstream &ofstr, const vector<vector<double>> &X) {
    size_t sizeX = X.size(), dimensionX = X[0].size();
    for (int i = 0; i < sizeX; i++) {
        for (int j = 0; j < dimensionX - 1; j++) {
            ofstr << X[i][j] << " ";
        }
        ofstr << X[i][dimensionX - 1] << "\n";
    }
    ofstr << "\n" << endl;
}

void addPointsGnuplot(ofstream &ofstr, const vector<TrialConstrained> &trials) {
    size_t numberPoints = trials.size();
    for (int i = 0; i < numberPoints; i++) {
        ofstr << trials[i].x << " " << trials[i].z << "\n";
    }
    ofstr << "\n" << endl;
}

void addPointsGnuplot(ofstream &ofstr, const vector<vector<double>> &X, const vector<TrialConstrained> &trials) {
    size_t sizeX = X.size(), dimensionX = X[0].size();
    for (int i = 0; i < sizeX; i++) {
        for (int j = 0; j < dimensionX; j++) {
            ofstr << X[i][j] << " ";
        }
        ofstr << trials[i].z << "\n";
    }
    ofstr << "\n" << endl;
} */
