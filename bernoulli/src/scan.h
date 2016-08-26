#ifndef SCANH
#define SCANH

void getCCCount(double * x, double * y, int * nCass, int * nCons, int locCount, double wSize, int wCount, int * casInW, int * conInW);

void loglikelihood(double * ll, int * casInW, int * conInW, int totalWindow, int casCount, int conCount, int highLow);

void findTopNCluster(double * x, double * y, int locCount, double * ll, double wSize, int wCount, int * center, int * radius,  double * cLL, int nClusters);

#endif
