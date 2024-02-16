#include <gnuplot/OutputFile.h>

void OutputFile::initArray(const std::string &name, size_t numberValues) {
    std::stringstream strArray;
    strArray << "array " << name << "[" << numberValues << "]\n";
    stream << strArray.str();
}

void OutputFile::addPoint(double x, double f, bool space) {
    std::stringstream strPoint;
    strPoint << x << " " << f << "\n";
    if (space) strPoint << "\n\n";
    stream << strPoint.str();
}

void OutputFile::addPoint(size_t a, double b, bool space) {
    std::stringstream strPoint;
    strPoint << a << " " << b << "\n";
    if (space) strPoint << "\n\n";
    stream << strPoint.str();
}

void OutputFile::addPoint(const std::vector<double> &x, double f, bool space) {
    std::stringstream strPoint;
    size_t dimensionalX = x.size();
    for (size_t i = 0; i < dimensionalX; ++i) {
        strPoint << x[i] << " ";
    }
    strPoint << f << "\n";
    if (space) strPoint << "\n\n";
    stream << strPoint.str();
}

void OutputFile::addPoints(const std::vector<Trial> &trials, bool space) {
    std::stringstream strPoints;
    size_t numberPoints = trials.size();
    for (size_t i = 0; i < numberPoints; ++i) {
        strPoints << trials[i].x << " " << trials[i].z << "\n";
    }
    if (space) strPoints << "\n\n";
    stream << strPoints.str();
}

void OutputFile::addPoints(const std::vector<opt::IndexTrial> &trials, bool space) {
    std::stringstream strPoints;
    size_t numberPoints = trials.size();
    for (size_t i = 0; i < numberPoints; ++i) {
        strPoints << trials[i].x << " " << trials[i].z << trials[i].nu << "\n";
    }
    if (space) strPoints << "\n\n";
    stream << strPoints.str();
}

void OutputFile::addPoints(const std::vector<double> &x, double f, bool space) {
    std::stringstream strPoints;
    size_t numberPoints = x.size();
    for (size_t i = 0; i < numberPoints; ++i) {
        strPoints << x[i] << " " << f << "\n";
    }
    if (space) strPoints << "\n\n";
    stream << strPoints.str();
}

void OutputFile::addPoints(const std::vector<std::vector<double>> &x, bool space) {
    std::stringstream strPoints;
    size_t numberPoints = x.size();
    if (numberPoints != 0) {
        size_t dimensionPoint = x[0].size();
        for (size_t i = 0; i < numberPoints; ++i) {
            for (size_t j = 0; j < dimensionPoint; ++j) {
                strPoints << x[i][j] << " ";
            }
            strPoints << "\n";
        }
        if (space) strPoints << "\n\n";
        stream << strPoints.str();
    }
}

/*
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
