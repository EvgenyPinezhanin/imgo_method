#ifndef SCRIPT_H_
#define SCRIPT_H_

#include <string>
#include <vector>

class Script {
private:
    std::string name, args;
    int error;

public:
    Script() : name(), args(), error(0) {};
    Script(const std::string &_name)
        : name(_name), args(), error(0) {};

    void setName(const std::string &_name) { name = _name; };

    int isError() const { return error; };

    template <typename T>
    void addArg(const T &arg) { args += std::to_string(arg) + " "; }

    void addArg(const std::string &arg) { args += "\"" + arg + "\" "; }

    template <typename T>
    void addArgs(const std::vector<T> &_args) {
        for (T arg : _args) addArg(arg);
    }

    void start();
};

#endif // SCRIPT_H_
