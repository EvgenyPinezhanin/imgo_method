#include <gnuplot/start_script.h>

using std::to_string;

void startScript(string nameScript) {
    startScript(nameScript, vector<int>{});
}

void startScript(string nameScript, int arg) {
    startScript(nameScript, vector<int>{ arg });
}

void startScript(string nameScript, const vector<int> &args) {
    string inputString;
    errorScript = 0;

#if defined(__linux__)
    setenv("QT_QPA_PLATFORM", "xcb", false);
    inputString = "chmod +x " + nameScript;
    if (system(inputString.c_str())) errorScript = 1;
#endif

    inputString = "gnuplot -c " + nameScript;
    for (int i = 0; i < args.size(); i++) {
        inputString += " " + to_string(args[i]);
    }
    if (system(inputString.c_str())) errorScript = 2;
}
