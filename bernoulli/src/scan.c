#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

void getCCCount(double * x, double * y, int * nCass, int * nCons, int locCount, double wSize, int wCount, int * casInW, int * conInW) {
	double distance;
	int minWindow;

#pragma omp parallel for private(distance, minWindow)
	for(int i = 0; i < locCount; i++) {
		for(int j = 0; j < locCount; j++) {
			distance = sqrt((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]));
			minWindow = (int)(ceil(distance / wSize));
			if(minWindow > 0)
				minWindow --;
			for(int k = minWindow; k < wCount; k++) {
				casInW[i * wCount + k] += nCass[j];
				conInW[i * wCount + k] += nCons[j];
			}
		}
	}

	return;
}

void loglikelihood(double * ll, int * casInW, int * conInW, int totalWindow, int casCount, int conCount, int highLow) {
	double cas, con, tot;
	double llTemp;
	int totCount = casCount + conCount;
	bool highCluster = true;
	bool lowCluster = true;
	if(highLow == 1)
		lowCluster = false;
	else if(highLow == -1)
		highCluster = false;

#pragma omp parallel for private(cas, con, tot, llTemp)
	for(int i = 0; i < totalWindow; i++) {
		cas = casInW[i];
		con = conInW[i];
		tot = cas + con;

		if(cas == -1) {
			ll[i] = 1;
		}
		else if(cas * conCount > con * casCount) { //High cluster of cases
			if(highCluster) {
				llTemp = cas * log(cas/tot);
				if(con > 0)
					llTemp += con * log(con/tot);
				if(casCount > cas)
					llTemp += (casCount - cas) * log((casCount - cas)/(totCount - tot));
				if(conCount > con)
					llTemp += (conCount - con) * log((conCount - con)/(totCount - tot));
				ll[i] = llTemp;
			}
			else
				ll[i] = 1;
		}
		else { //Low cluster of cases
			if(lowCluster) {
				llTemp = con * log(con/tot);
				if(cas > 0)
					llTemp += cas * log(cas/tot);
				if(casCount > cas)
					llTemp += (casCount - cas) * log((casCount - cas)/(totCount - tot));
				if(conCount > con)
					llTemp += (conCount - con) * log((conCount - con)/(totCount - tot));
				ll[i] = llTemp;
			}
			else
				ll[i] = 1;
		}
	}

	return;	
}

void findTopNCluster(double * x, double * y, int locCount, double * ll, double wSize, int wCount, int * center, int * radius,  double * cLL, int nClusters) {
	if(nClusters < 1)
		return;
	
	int aCenter = -1;
	int aRadius = -1;
	
	for(int i = 0; i < locCount; i++) {
		for(int j = 0; j < wCount; j++) {
			if(ll[i * wCount + j] < 0) {
				if(aCenter < 0) {
					aCenter = i;
					aRadius = j;
				}
				else if(ll[i * wCount + j] > ll[aCenter * wCount + aRadius]) {
					aCenter = i;
					aRadius = j;
				}
			}
		}
	}
	
	center[0] = aCenter;
	radius[0] = aRadius;
	cLL[0] = ll[aCenter * wCount + aRadius];

	double lastX, lastY, lastRad;
	lastX = x[aCenter];
	lastY = y[aCenter];
	lastRad = (aRadius + 1) * wSize;

	double distance;
	int maxWindow;

	for(int c = 1; c < nClusters; c ++) {
		//Remove intersecting clusters
		for(int i = 0; i < locCount; i++) {
			distance = sqrt((x[i] - lastX) * (x[i] - lastX) + (y[i] - lastY) * (y[i] - lastY)) - lastRad;
			maxWindow = ceil(distance / wSize) - 1;
			if(maxWindow < 0)
				maxWindow = 0;
			for(int j = maxWindow; j < wCount; j++)
				ll[i * wCount + j] = 1;
		}

		//Find secoundary clusters
		aCenter = -1;
		aRadius = -1;

		for(int i = 0; i < locCount; i++) {
			for(int j = 0; j < wCount; j++) {
				if(ll[i * wCount + j] < 0) {
					if(aCenter < 0) {
						aCenter = i;
						aRadius = j;
					}
					else if(ll[i * wCount + j] > ll[aCenter * wCount + aRadius]) {
						aCenter = i;
						aRadius = j;
					}
				}
			}
		}
		center[c] = aCenter;
		radius[c] = aRadius;
		if(aCenter != -1)
			cLL[c] = ll[aCenter * wCount + aRadius];
		else
			break;
		lastX = x[aCenter];
		lastY = y[aCenter];
		lastRad = (aRadius + 1) * wSize;						
	}
	
	return;
}
