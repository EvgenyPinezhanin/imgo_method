#ifndef TRIALS_FILE_H_
#define TRIALS_FILE_H_

#include <fstream>
#include <string>
#include <vector>

#include <trials/Trial.h>

using std::ofstream;
using std::string;
using std::vector;

class TrialsFile {
private:
    ofstream file;
    bool open;

public:
    TrialsFile() : file(), open(false) {};
    TrialsFile(const string &fileName) : file(fileName), open(file.is_open()) {};

    void newFile(const string &fileName);
    void closeFile();

    bool isOpen() const { return open; };

    void addPoint(double x, double f, bool space = true);

    // void addPoint(ofstream &ofstr, Trial trial, bool space = true);
    // void addPointGnuplot(ofstream &ofstr, const vector<double> &x, double f);
    // void addPointGnuplot(ofstream &ofstr, const vector<double> &x);

    void addPoints(const vector<Trial> &trials, bool space = true);

    // void addPointsGnuplot(ofstream &ofstr, const vector<vector<double>> &X, const vector<double> &f);
    // void addPointsGnuplot(ofstream &ofstr, const vector<vector<double>> &X);
    // void addPointsGnuplot(ofstream &ofstr, const vector<double> &X, const vector<double> &f);
    // void addPointsGnuplot(ofstream &ofstr, const vector<TrialConstrained> &trials);
    // void addPointsGnuplot(ofstream &ofstr, const vector<vector<double>> &X, const vector<TrialConstrained> &trials);
};

#endif // TRIALS_FILE_H_
