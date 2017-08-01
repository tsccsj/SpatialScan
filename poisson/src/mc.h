#ifndef MCH
#define MCH

int * monteCarlo(double * x, double * y, double * intensity, double * intenInW, int locCount, int casCount, double wSize, int wCount, bool highLow, double * clusterLL, int nClusters, int nSim);
int * monteCarloOld(double * x, double * y, double * intensity, int locCount, int casCount, int * clusterCase, int * centerID, double * cRadius, bool * highCluster, int nClusters, int nSim);

#endif
