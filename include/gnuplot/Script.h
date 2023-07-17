#ifndef SCRIPT_H_
#define SCRIPT_H_

#include <string>
#include <vector>

using std::string;
using std::to_string;
using std::vector;

class Script {
private:
    string name, args;
    int error;

public:
    Script() : name(), error(0) {};
    Script(const string &_name) : name(_name), error(0) {};

    void setName(const string &_name) { name = _name; };

    int isError() const { return error; };

    template <typename T>
    void addArg(const T &arg) { args += to_string(arg) + " "; }

    void addArg(const string &arg) { args += "\"" + arg + "\" "; }

    template <typename T>
    void addArgs(const vector<T> &_args) { for (T arg : _args) addArg(arg); }

    void start();
};

#endif // SCRIPT_H_
