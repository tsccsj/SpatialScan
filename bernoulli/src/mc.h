#ifndef MCH
#define MCH

int * monteCarlo(double * x, double * y, int * locEnding, int locCount, int casCount, int allCount, double wSize, int wCount, int highLow, double * clusterLL, int nClusters, int nSim);

int * monteCarloOld(double * x, double * y, int * locEnding, int locCount, int casCount, int allCount, int * clusterCase, int * centerID, double * cRadius, bool * highCluster, int nClusters, int nSim);

#endif
