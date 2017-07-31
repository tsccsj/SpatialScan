#ifndef MCH
#define MCH

int * monteCarlo(double * x, double * y, int * locEnding, int locCount, int casCount, int allCount, int wSize, int wCount, bool highLow, double * clusterLL, bool * highCluster, int nClusters, int nSim);

int * monteCarloOld(double * x, double * y, int * locEnding, int locCount, int casCount, int allCount, int * clusterCase, int * centerID, double * cRadius, bool * highCluster, int nClusters, int nSim);

#endif
