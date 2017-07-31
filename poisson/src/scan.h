#ifndef SCANH
#define SCANH

void getWindowCandI(double * x, double * y, int * nCass, double * intensity, int locCount, double wSize, int wCount, int * casInW, double * intenInW);

void getWindowCOnly(double * x, double * y, int * nCass, int locCount, double wSize, int wCount, int * casInW);

void loglikelihood(double * ll, int * casInW, double * intenInW, int totalWindow, int casCount, int highLow);

void findTopNCluster(double * x, double * y, int locCount, double * ll, double wSize, int wCount, int * center, int * radius,  double * cLL, int nClusters);

#endif
