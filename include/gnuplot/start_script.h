#ifndef START_SCRIPT_H_
#define START_SCRIPT_H_

#include <string>
#include <vector>

using std::string;
using std::vector;

static int errorScript;

void startScript(string nameScript);
void startScript(string nameScript, int arg);
void startScript(string nameScript, const vector<int> &args);

#endif // START_SCRIPT_H_
