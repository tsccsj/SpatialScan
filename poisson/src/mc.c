#include <stdio.h>
#include <stdlib.h>
#include <random>
#include <omp.h>

using namespace std;

void simulateCases(double * preInten, int * simCases, int locCount, int casCount) {
	
	static std::random_device rd;
	static std::mt19937 rng(rd());
	static std::uniform_real_distribution<double> uni(0, preInten[locCount - 1]); //[a, b)

	for(int i = 0; i < casCount; i++) {
		simCases[i] = 0;
	}

	double randomNumber;
	for(int i = 0; i < casCount; i++) {
		randomNumber = uni(rng);
		for(int j = 0; j < locCount; j++) {
			if(randomNumber < preInten[j]) {
				simCases[j] ++;
				break;
			}
		}	
	}

	
}

int * monteCarlo(double * x, double * y, double * intensity, int locCount, int casCount, int * clusterCase, int * centerID, double * cRadius, bool * highCluster, int nClusters, int nSim) {
	int * nExtreme;
	
	if(NULL == (nExtreme = (int *) malloc (nClusters * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}


	double * preInten;
	int * simCass;

	if(NULL == (preInten = (double *) malloc (locCount * sizeof(double)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	
	preInten[0] = intensity[0];

	for(int i = 1; i < locCount; i++) {
		preInten[i] = preInten[i-1] + intensity[i];
	}

	if(NULL == (simCass = (int *) malloc (locCount * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	for(int i = 0; i < nSim; i++) {

		simulateCases(preInten, simCass, locCount, casCount);
		
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


	free(preInten);
	free(simCass);

	return nExtreme;	

}
