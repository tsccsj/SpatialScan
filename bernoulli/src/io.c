/**
 * io.c
 * Authors: Yizhao Gao <yizhaotsccsj@gmail.com>
 * Date: {08/01/2017}
 */

#include <stdio.h>
#include <stdlib.h>

int getNumPoints(FILE * file) {
	rewind(file);
	int count = 0;

	double x, y;
	int nCas, nCon;

	while(EOF != fscanf(file, "%lf,%lf,%d,%d\n", &x, &y, &nCas, &nCon)) {
		count ++;
	}

	return count;
}

void readFile(FILE * file, double * x, double * y, int * nCass, int * nCons, int & casCount, int & conCount) {
	rewind(file);
	
	int locID = 0;

	casCount = 0;
	conCount = 0;

	while(EOF != fscanf(file, "%lf,%lf,%d,%d\n", x + locID, y + locID, nCass + locID, nCons + locID)) {
		casCount += nCass[locID];
		conCount += nCons[locID];
		locID ++;
	}

	return;
}
