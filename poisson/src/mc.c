/**
 * mc.c
 * Authors: Yizhao Gao <yizhaotsccsj@gmail.com>
 * Date: {08/01/2017}
 */

#include <stdio.h>
#include <stdlib.h>
#include <random>
#include <omp.h>
#include "scan.h"

using namespace std;

void simulateCases(double * intensity, int * simCases, int locCount, int casCount) {
	
	static std::random_device rd;
	static std::mt19937 rng(rd());
	static std::discrete_distribution<int> d (intensity, intensity + locCount);

	for(int i = 0; i < locCount; i++) {
		simCases[i] = 0;
	}

	for(int i = 0; i < casCount; i++) {
		simCases[d(rng)] ++;
	}

	
}

int * monteCarlo(double * x, double * y, double * intensity, double * intenInW, int locCount, int casCount, double wSize, int wCount, bool highLow, double * clusterLL, int nClusters, int nSim) {

	int * nExtreme;

	if(NULL == (nExtreme = (int *) malloc (nClusters * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	for(int i = 0; i < nClusters; i++)
		nExtreme[i] = 0;

	int * simCass;
	if(NULL == (simCass = (int *) malloc (locCount * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	int * simCasInW;
	double * simll;

	if(NULL == (simCasInW = (int *) malloc (locCount * wCount * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (simll = (double *) malloc (locCount * wCount * sizeof(double)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	double simMaxLL;

	for(int i = 0; i < nSim; i++) {
		simulateCases(intensity, simCass, locCount, casCount);

		getWindowCOnly(x, y, simCass, locCount, wSize, wCount, simCasInW);
		
		loglikelihood(simll, simCasInW, intenInW, locCount * wCount, casCount, highLow);
		
		simMaxLL = -9999;
		int k = 0;
		for(; k < locCount * wCount; k++) {
			if(simll[k] > 0) {
				simMaxLL = simll[k];
				k++;
				break;
			}
		}

		for(; k < locCount * wCount; k++) {
			if(simll[k] > 0 && simll[k] > simMaxLL) {
				simMaxLL = simll[k];
			}
		}
	
		if(simMaxLL > 0) {
			for(int j = 0; j < nClusters; j++) {
				if(simMaxLL > clusterLL[j]) {
					nExtreme[j] ++;
				}			
			}
		}
			
	}	

	free(simCasInW);
	free(simll);

	free(simCass);

	return nExtreme;
}


int * monteCarloOld(double * x, double * y, double * intensity, int locCount, int casCount, int * clusterCase, int * centerID, double * cRadius, bool * highCluster, int nClusters, int nSim) {
	int * nExtreme;
	
	if(NULL == (nExtreme = (int *) malloc (nClusters * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	for(int i = 0; i < nClusters; i++) {
		nExtreme[i] = 0;
	}

	int * simCass;

	if(NULL == (simCass = (int *) malloc (locCount * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	for(int i = 0; i < nSim; i++) {

		simulateCases(intensity, simCass, locCount, casCount);
		
#pragma omp parallel for
		for(int j = 0; j < nClusters; j++) {
			double xC = x[centerID[j]];
			double yC = y[centerID[j]];
			double rad2 = cRadius[j] * cRadius[j];
			int simCasInc = 0;
			for(int k = 0; k < locCount; k++) {
				if((x[k] - xC) * (x[k] - xC) + (y[k] - yC) * (y[k] - yC) <= rad2) {
					simCasInc += simCass[k];
				}	
			}
			if(highCluster[j] && simCasInc >= clusterCase[j])
				nExtreme[j] ++;
			else if(!highCluster[j] && simCasInc <= clusterCase[j])
				nExtreme[j] ++;
		}
	}


	free(simCass);

	return nExtreme;	

}
