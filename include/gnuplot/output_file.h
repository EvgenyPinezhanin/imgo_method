#ifndef OUTPUT_FILE_H_
#define OUTPUT_FILE_H_

#include <fstream>
#include <vector>
#include <string>
#include <sstream>

#include <trials/Trial.h>
#include <base_structures/trials/IndexTrial.h>

class output_file {
private:
    std::ofstream stream;

public:
    output_file() : stream() {};
    output_file(const std::string &name) : stream(name) {};

    void open(const std::string &name) { stream.open(name); };
    void close() { stream.close(); };

    bool is_open() const { return stream.is_open(); };

    template<typename T>
    void set_variable(const std::string &name, T value, bool is_string_expression = true) {
        std::stringstream str_variable;
        str_variable << name << (is_string_expression ? " = \"" : " = ") << value << (is_string_expression ? "\"\n" : "\n");
        stream << str_variable.str();
    }

    void init_array(const std::string &name, size_t number_values);

    template<typename T>
    void set_value_in_array(const std::string &name, size_t index, T value, bool is_string_expression = true) {
        std::stringstream str_value;
        str_value << name << "[" << index << (is_string_expression ? "] = \"" : "] = ") <<
                                    value << (is_string_expression ? "\"\n" : "\n");
        stream << str_value.str();
    }

    template<typename T>
    void set_values_in_array(const std::string &name, size_t start_index, std::vector<T> &values, bool is_string_expression = true) {
        size_t last_index = start_index + values.size();
        for (size_t i = start_index; i < last_index; i++) {
            set_value_in_array(name, i, values[i - start_index], is_string_expression);
        }
    }

    template <typename T>
    void init_array(const std::string &name, const std::vector<T> &values, bool is_string_expression = true) {
        size_t number_values = values.size();
        init_array(name, number_values);
        for (size_t i = 0; i < number_values; ++i) {
            set_value_in_array(name, i + 1, values[i], is_string_expression);
        }
    }

    void add_point(double x, double f, bool space = true);
    void add_point(const std::vector<double> &x, double f, bool space = true);

    // void addPoint(ofstream &ofstr, Trial trial, bool space = true);
    // void addPointGnuplot(ofstream &ofstr, const vector<double> &x, double f);

    void add_points(const std::vector<Trial> &trials, bool space = true);
    void add_points(const std::vector<opt::IndexTrial> &trials, bool space = true);
    void add_points(const std::vector<double> &x, double f, bool space = true);

    // void addPointsGnuplot(ofstream &ofstr, const vector<vector<double>> &X, const vector<double> &f);
    // void addPointsGnuplot(ofstream &ofstr, const vector<vector<double>> &X);
    // void addPointsGnuplot(ofstream &ofstr, const vector<double> &X, const vector<double> &f);
    // void addPointsGnuplot(ofstream &ofstr, const vector<TrialConstrained> &trials);
    // void addPointsGnuplot(ofstream &ofstr, const vector<vector<double>> &X, const vector<TrialConstrained> &trials);
};

#endif // OUTPUT_FILE_H_
