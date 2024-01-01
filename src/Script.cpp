#include <gnuplot/Script.h>

void Script::start() {
    std::string inputString;
    error = 0;

#if defined(__linux__)
    setenv("QT_QPA_PLATFORM", "xcb", false);
    inputString = "chmod +x " + name;
    if (system(inputString.c_str())) error = 1;
#endif

    inputString = "gnuplot -c " + name + " " + args;
    if (system(inputString.c_str())) error = 2;
}
