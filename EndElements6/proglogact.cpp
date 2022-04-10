#include <string>
#include <fstream>
#include "proglogact.h"

using namespace std;

void writeErr(string error)
{
	ofstream errFile;
	errFile.open("Error.txt");
	errFile << error << endl;
	errFile.close();
}