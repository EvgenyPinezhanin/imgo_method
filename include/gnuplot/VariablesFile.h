#ifndef VARIABLES_FILE_H_
#define VARIABLES_FILE_H_

#include <fstream>
#include <string>
#include <vector>

using std::ofstream;
using std::vector;
using std::string;
using std::endl;

class VariablesFile {
private:
    ofstream file;
    bool open;

public:
    VariablesFile() : file(), open(false) {};
    VariablesFile(const string &fileName) : file(fileName), open(file.is_open()) {};

    void newFile(const string &fileName);
    void closeFile();

    bool isOpen() const { return open; };

    template<typename T>
    void setVariable(string nameVariable, T value, bool stringExpression = true) {
        file << nameVariable << (stringExpression ? " = \"" : " = ") << value << (stringExpression ? "\"" : "") << endl;
    }

    void initArray(string nameArray, int numberElements);

    template<typename T>
    void setValueInArray(string nameArray, int index, T value, bool stringExpression = true) {
        file << nameArray << "[" << index << (stringExpression ? "] = \"" : "] = ") <<
                                    value << (stringExpression ? "\"" : "") << endl;
    }

    template <typename T>
    void initArray(string nameArray, vector<T> elements, bool stringExpression = true) {
        int numberElements = (int)elements.size();
        initArray(nameArray, numberElements);
        for (int i = 0; i < numberElements; i++) {
            setValueInArray(nameArray, i + 1, elements[i], stringExpression);
        }
    }
};

#endif // VARIABLES_FILE_H_
