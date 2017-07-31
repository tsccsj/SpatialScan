#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

void getWindowCandI(double * x, double * y, int * nCass, double * intensity, int locCount, double wSize, int wCount, int * casInW, double * intenInW) {
	double distance;
	int minWindow;

	for(int i = 0; i < locCount * wCount; i++) {
		casInW[i] = 0;
		intenInW[i] = 0.0;
	}

#pragma omp parallel for private(distance, minWindow)
	for(int i = 0; i < locCount; i++) {
		for(int j = 0; j < locCount; j++) {
			distance = sqrt((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]));
			minWindow = (int)(ceil(distance / wSize));
			if(minWindow > 0)
				minWindow --;
			for(int k = minWindow; k < wCount; k++) {
				casInW[i * wCount + k] += nCass[j];
				intenInW[i * wCount + k] += intensity[j];
			}
		}
	}

	return;
}

void getWindowCOnly(double * x, double * y, int * nCass, int locCount, double wSize, int wCount, int * casInW) {
	double distance;
	int minWindow;
	
	for(int i = 0; i < locCount * wCount; i++) {
		casInW[i] = 0;
	}

#pragma omp parallel for private(distance, minWindow)
	for(int i = 0; i < locCount; i++) {
		for(int j = 0; j < locCount; j++) {
			distance = sqrt((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]));
			minWindow = (int)(ceil(distance / wSize));
			if(minWindow > 0)
				minWindow --;
			for(int k = minWindow; k < wCount; k++) {
				casInW[i * wCount + k] += nCass[j];
			}
		}
	}

	return;
}

void loglikelihood(double * ll, int * casInW, double * intenInW, int totalWindow, int casCount, int highLow) {
	double cas, inten;
	double llTemp;
	bool highCluster = true;
	bool lowCluster = true;
	if(highLow == 1)
		lowCluster = false;
	else if(highLow == -1)
		highCluster = false;

#pragma omp parallel for private(cas, inten, llTemp)
	for(int i = 0; i < totalWindow; i++) {
		cas = casInW[i];
		inten = intenInW[i];

		if(cas == -1) {
			ll[i] = -9999;
		}
		else if(cas > inten) { //High cluster of cases
			if(highCluster) {
				llTemp = cas * log(cas/inten);
				if(cas < casCount)
					llTemp += (casCount - cas) * log((casCount - cas)/(casCount - inten));
				ll[i] = llTemp;
			}
			else
				ll[i] = -9999;
		}
		else if(cas < inten) { //Low cluster of cases
			if(lowCluster) {
				llTemp = (casCount - cas) * log((casCount - cas)/(casCount - inten));
				if(cas > 0)
					llTemp += cas * log(cas/inten);
				ll[i] = llTemp;
			}
			else
				ll[i] = -9999;
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
			if(ll[i * wCount + j] > -9990) {
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
				ll[i * wCount + j] = -9999;
		}

		//Find secoundary clusters
		aCenter = -1;
		aRadius = -1;

		for(int i = 0; i < locCount; i++) {
			for(int j = 0; j < wCount; j++) {
				if(ll[i * wCount + j] > -9990) {
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
