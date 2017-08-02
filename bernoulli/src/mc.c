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

void randomLabel(int * indAll, int casCount, int allCount) {
	static std::random_device rd;
	static std::mt19937 rng(rd());
	static std::uniform_int_distribution<int> uni(0, allCount - 1);

	int casID;
	for(int i = 0; i < allCount; i++)
		indAll[i] = 0;

	for(int i = 0; i < casCount; i++) {
		casID = uni(rng);
		while(indAll[casID] == 1)
			casID = uni(rng);
		indAll[casID] = 1;
	}

	return;
}

int * monteCarlo(double * x, double * y, int * locEnding, int locCount, int casCount, int allCount, double wSize, int wCount, bool highLow, double * clusterLL, int nClusters, int nSim) {

	int * nExtreme;
	
	if(NULL == (nExtreme = (int *) malloc (nClusters * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	for(int i = 0; i < nClusters; i++)
		nExtreme[i] = 0;

	int * indAll;
	int * simCass;
	int * simCons;

	if(NULL == (indAll = (int *) malloc (allCount * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (simCass = (int *) malloc (locCount * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (simCons = (int *) malloc (locCount * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	int indID, simCas, simCon;

	int * simCasInW;
	int * simConInW;
	double * simll;
	
	if(NULL == (simCasInW = (int *) malloc (locCount * wCount * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (simConInW = (int *) malloc (locCount * wCount * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (simll = (double *) malloc (locCount * wCount * sizeof(double)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	double simMaxLL;

	for(int i = 0; i < nSim; i++) {
		randomLabel(indAll, casCount, allCount);

		indID = 0;
		for(int j = 0; j < locCount; j++) {
			simCas = 0;
			simCon = 0;
			for(; indID < locEnding[j]; indID ++) {
				if(indAll[indID] == 1) {
					simCas ++;
				}
				else {
					simCon ++;
				}
			}
			simCass[j] = simCas;
			simCons[j] = simCon;
		}

		getCCCount(x, y, simCass, simCons, locCount, wSize, wCount, simCasInW, simConInW);	
		loglikelihood(simll, simCasInW, simConInW, locCount * wCount, casCount, allCount - casCount, highLow);
		
		simMaxLL = 1;
		int k = 0;
		for(; k < locCount * wCount; k++) {
			if(simll[k] < 0) {
				simMaxLL = simll[k];
				k++;
				break;
			}
		}

		for(; k < locCount * wCount; k++) {
			if(simll[k] < 0 && simll[k] > simMaxLL) {
				simMaxLL = simll[k];
			}
		}
	
		if(simMaxLL < 0) {
			for(int j = 0; j < nClusters; j++) {
				if(simMaxLL > clusterLL[j]) {
					nExtreme[j] ++;
				}			
			}
		}
		
	}


	free(simCasInW);
	free(simConInW);
	free(simll);

	free(indAll);
	free(simCass);
	free(simCons);
	
	return nExtreme;
}


int * monteCarloOld(double * x, double * y, int * locEnding, int locCount, int casCount, int allCount, int * clusterCase, int * centerID, double * cRadius, bool * highCluster, int nClusters, int nSim) {
	int * nExtreme;

	if(NULL == (nExtreme = (int *) malloc (nClusters * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	for(int i = 0; i < nClusters; i++)
		nExtreme[i] = 0;

	int * indAll;
	int * simCass;
	if(NULL == (indAll = (int *) malloc (allCount * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (simCass = (int *) malloc (locCount * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	int indID, simCas;

	for(int i = 0; i < nSim; i++) {
		randomLabel(indAll, casCount, allCount);

		indID = 0;
		for(int j = 0; j < locCount; j++) {
			simCas = 0;
			for(; indID < locEnding[j]; indID ++) {
				if(indAll[indID] == 1)
					simCas ++;
			}
			simCass[j] = simCas;
		}

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

	free(indAll);
	free(simCass);

	return nExtreme;

}


