#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mc.h"
#include "scan.h"
#include "io.h"


int main(int argc, char ** argv) {

	if(argc != 7) {
		printf("Incorrect number of parameters: SpatialScanPoisson [inputFile] [windowInc] [windowCount] [#ofClusters] [#ofMonteCarlo] [HighOrLowIndicator]\n");
		printf("[HighOrLowIndicator]\n\t1: High Only\n\t-1: Low Only\n\t0: Both\n");
		exit(1);
	}

	double wSize = atof(argv[2]);
	int wCount = atoi(argv[3]);
	int nClusters = atoi(argv[4]);
	int nSim = atoi(argv[5]);
	int highLow = atoi(argv[6]);

	double * x;
	double * y;
	int * nCass;
	double * intensity;
	
	int casCount;
	double totalInten;
	int locCount;
	
	FILE * file;
	
	if(NULL == (file = fopen(argv[1], "r"))) {
		printf("ERROR: Can't open the input file.\n");
		exit(1);
	}
	

	locCount = getNumPoints(file);

	if(NULL == (x = (double *) malloc (locCount * sizeof(double)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (y = (double *) malloc (locCount * sizeof(double)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (nCass = (int *) malloc (locCount * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (intensity = (double *) malloc (locCount * sizeof(double)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	readFile(file, x, y, nCass, intensity, casCount, totalInten);
	
	printf("There is %d locations\n", locCount);
	printf("Total CASE count: \t%d\nTotal raw intensity:\t%lf\n", casCount, totalInten);

	for(int i = 0; i < locCount; i++) {
		intensity[i] = intensity[i] * casCount / totalInten;
	}

	fclose(file);
	printf("Finish reading input files.\n");
	
	int * casInW;
	double * intenInW;

	if(NULL == (casInW = (int *) malloc (locCount * wCount * sizeof(int)))) {	
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	if(NULL == (intenInW = (double *) malloc (locCount * wCount * sizeof(double)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	for(int i = 0; i < locCount * wCount; i++) {
		casInW[i] = 0;
		intenInW[i] = 0.0;
	}

	getWindowCandI(x, y, nCass, intensity, locCount, wSize, wCount, casInW, intenInW);

	double * ll;
	if(NULL == (ll = (double *) malloc (locCount * wCount * sizeof(double)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	loglikelihood(ll, casInW, intenInW, locCount * wCount, casCount, highLow);

	int * center;
	int * radius;
	double * cLL;


	if(NULL == (center = (int *) malloc (nClusters * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	if(NULL == (radius = (int *) malloc (nClusters * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	
	if(NULL == (cLL = (double *) malloc (nClusters * sizeof(double)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	findTopNCluster(x, y, locCount, ll, wSize, wCount, center, radius,  cLL, nClusters);
	
	int * clusterCas;
	double * clusterInten;
	double * cRadius;
	bool * highCluster;

	if(NULL == (clusterCas = (int *) malloc (nClusters * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (clusterInten = (double *) malloc (nClusters * sizeof(double)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (cRadius = (double *) malloc (nClusters * sizeof(double)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (highCluster = (bool *) malloc (nClusters * sizeof(bool)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	
	int aCenter;
	int aRadius;

	for(int i = 0; i < nClusters; i ++) {
		aCenter = center[i];
		aRadius = radius[i];
		if(aCenter == -1) {
			nClusters = i;
			break;
		}
		clusterCas[i] = casInW[aCenter * wCount + aRadius];
		clusterInten[i] = intenInW[aCenter * wCount + aRadius];
		cRadius[i] = wSize * (aRadius + 1);

		if(clusterCas[i] > clusterInten[i]) 	//High clusters
			highCluster[i] = true;
		else
			highCluster[i] = false;
	}

	
//Here is the Monte Carlo Simulation
	int * nExtreme;

	if(nSim > 0) {
		nExtreme = monteCarlo(x, y, intensity, locCount, casCount, clusterCas, center, cRadius, highCluster, nClusters, nSim);
	}


	printf("############### Cluster Info ###############\n");
	printf("ID,HL,X,Y,Radius,#Cas,Exp#Cas");
	if(nSim > 0)
		printf(",LL,P\n");
	else
		printf(",LL\n");

	for(int i = 0; i < nClusters; i ++) {
		aCenter = center[i];
		aRadius = radius[i];

		printf("%d", i);
		if(highCluster[i]) 		//High clusters
			printf(",H");
		else				//Low clusters
			printf(",L");
		printf(",%lf,%lf,%lf", x[aCenter], y[aCenter], cRadius[i]);
		printf(",%d,%lf", clusterCas[i], clusterInten[i]);
		if(nSim > 0)
			printf(",%lf,%lf\n", cLL[i], (double)(nExtreme[i] + 1) / (nSim + 1));
		else
			printf(",%lf\n", cLL[i]);
	}

	printf("############ Cluster Membership ############\n");

	bool inCluster;
	double distance;		
	double radiusValue;

	//Find the cluster belonging of each location
	for(int i = 0; i < locCount; i++) {
		inCluster = false;
		for(int j = 0; j < nClusters; j++) {
			aCenter = center[j];
			radiusValue = cRadius[j];
			distance = sqrt((x[i] - x[aCenter]) * (x[i] - x[aCenter]) + (y[i] - y[aCenter]) * (y[i] - y[aCenter]));
			if(distance <= radiusValue) {
				printf("%lf,%lf,%d,%lf,%d\n", x[i], y[i], nCass[i], intensity[i], j);
				inCluster = true;
				break;
			}
		}
		if(!inCluster)
			printf("%lf,%lf,%d,%lf,-1\n", x[i], y[i], nCass[i], intensity[i]);
	}

	free(clusterCas);
	free(clusterInten);
	free(cRadius);
	free(highCluster);

	free(center);
	free(radius);
	free(cLL);

	free(ll);

	free(casInW);
	free(intenInW);

	free(x);
	free(y);
	free(nCass);
	free(intensity);

	if(nSim > 0) {
		free(nExtreme);
	}

	return 0;
}
