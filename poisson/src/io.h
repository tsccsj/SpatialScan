#ifndef IOH
#define IOH

int getNumPoints(FILE * file);

void readFile(FILE * file, double * x, double * y, int * nCass, double * intensity, int & casCount, double & totalInten);

#endif
