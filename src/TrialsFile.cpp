#include <gnuplot/TrialsFile.h>

using std::endl;

void TrialsFile::newFile(const string &fileName) {
    if (open) file.close();
    file.open(fileName);
    open = file.is_open();
}

void TrialsFile::closeFile() {
    file.close();
    open = file.is_open();
}

void TrialsFile::addPoint(double x, double f, bool space) {
    file << x << " " << f << "\n";
    if (space) file << "\n\n";
}

void TrialsFile::addPoints(const vector<Trial> &trials, bool space) {
    int numberPoints = (int)trials.size();
    for (int i = 0; i < numberPoints; i++) {
        file << trials[i].x << " " << trials[i].z << "\n";
    }
    if (space) file << "\n\n";
}

void TrialsFile::addPoints(const vector<double> &x, double f, bool space) {
    int numberPoints = (int)x.size();
    for (int i = 0; i < numberPoints; i++) {
        file << x[i] << " " << f << "\n";
    }
    if (space) file << "\n\n";
}

/* 

void addPointGnuplot(ofstream &ofstr, const vector<double> &X, double f) {
    size_t dimensionX = X.size();
    for (int i = 0; i < dimensionX; i++) {
        ofstr << X[i] << " ";
    }
    ofstr << f << "\n";
    ofstr << "\n" << endl;
}

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