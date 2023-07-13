#include <gnuplot/VariablesFile.h>

void VariablesFile::newFile(const string &fileName) {
    if (open) file.close();
    file.open(fileName);
    open = file.is_open();
}

void VariablesFile::closeFile() {
    file.close();
    open = file.is_open();
}

void VariablesFile::initArray(string nameArray, int numberElements) {
    file << "array " << nameArray << "[" << numberElements << "]\n";
}
