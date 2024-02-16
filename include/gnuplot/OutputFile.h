#ifndef OUTPUT_FILE_H_
#define OUTPUT_FILE_H_

#include <fstream>
#include <vector>
#include <string>
#include <sstream>

#include <trials/Trial.h>
#include <base_structures/trials/IndexTrial.h>

class OutputFile {
private:
    std::ofstream stream;

public:
    OutputFile() : stream() {};
    OutputFile(const std::string &name) : stream(name) {};

    void open(const std::string &name) { stream.open(name); };
    void close() { stream.close(); };

    bool isOpen() const { return stream.is_open(); };

    template<typename T>
    void setVariable(const std::string &name, T value, bool isString = true) {
        std::stringstream strVariable;
        strVariable << name << (isString ? " = \"" : " = ") << value << (isString ? "\"\n" : "\n");
        stream << strVariable.str();
    }

    void initArray(const std::string &name, size_t numberValues);

    template<typename T>
    void setValueInArray(const std::string &name, size_t index, T value, bool isString = true) {
        std::stringstream strValue;
        strValue << name << "[" << index << (isString ? "] = \"" : "] = ") <<
                                   value << (isString ? "\"\n" : "\n");
        stream << strValue.str();
    }

    template<typename T>
    void setValuesInArray(const std::string &name, size_t startIndex, std::vector<T> &values, bool isString = true) {
        size_t lastIndex = startIndex + values.size();
        for (size_t i = startIndex; i < lastIndex; ++i) {
            setValueInArray(name, i, values[i - startIndex], isString);
        }
    }

    template <typename T>
    void initArray(const std::string &name, const std::vector<T> &values, bool isString = true) {
        size_t numberValues = values.size();
        initArray(name, numberValues);
        for (size_t i = 0; i < numberValues; ++i) {
            setValueInArray(name, i + 1, values[i], isString);
        }
    }

    void addPoint(double x, double f, bool space = true);
    void addPoint(size_t a, double b, bool space = true);
    void addPoint(const std::vector<double> &x, double f, bool space = true);

    // void addPoint(ofstream &ofstr, Trial trial, bool space = true);
    // void addPointGnuplot(ofstream &ofstr, const vector<double> &x, double f);

    void addPoints(const std::vector<Trial> &trials, bool space = true);
    void addPoints(const std::vector<opt::IndexTrial> &trials, bool space = true);
    void addPoints(const std::vector<double> &x, double f, bool space = true);
    void addPoints(const std::vector<std::vector<double>> &x, bool space = true);

    // void addPointsGnuplot(ofstream &ofstr, const vector<vector<double>> &X, const vector<double> &f);
    // void addPointsGnuplot(ofstream &ofstr, const vector<double> &X, const vector<double> &f);
    // void addPointsGnuplot(ofstream &ofstr, const vector<TrialConstrained> &trials);
    // void addPointsGnuplot(ofstream &ofstr, const vector<vector<double>> &X, const vector<TrialConstrained> &trials);
};

#endif // OUTPUT_FILE_H_
