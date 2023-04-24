#ifndef TEXTPARSER_H
#define TEXTPARSER_H

#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <ctime>
#include <fstream>
#include <iomanip>

using namespace std;

int parse(FILE *stream, vector<string> &taxa, vector<string> &sequences);
int writeDmat(double **DMat, vector<string> &taxa);

#endif
