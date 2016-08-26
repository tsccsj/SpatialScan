#ifndef IOH
#define IOH

int getNumPoints(FILE * file);

void readFile(FILE * file, double * x, double * y, int * nCass, int * nCons, int & casCount, int & conCount);

#endif
