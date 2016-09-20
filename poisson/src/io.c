#include <stdio.h>
#include <stdlib.h>

int getNumPoints(FILE * file) {
	rewind(file);
	int count = 0;

	double x, y;
	int nCas;
	double inten;

	while(EOF != fscanf(file, "%lf,%lf,%d,%lf\n", &x, &y, &nCas, &inten)) {
		count ++;
	}

	return count;
}

void readFile(FILE * file, double * x, double * y, int * nCass, double * intensity, int & casCount, double & totalInten) {
	rewind(file);
	
	int locID = 0;

	casCount = 0;
	totalInten = 0;

	while(EOF != fscanf(file, "%lf,%lf,%d,%lf\n", x + locID, y + locID, nCass + locID, intensity + locID)) {
		casCount += nCass[locID];
		totalInten += intensity[locID];
		locID ++;
	}

	return;
}
