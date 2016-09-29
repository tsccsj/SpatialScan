# Spatial scan with Poisson Model
Used for finding spatial clusters of events over a background intensity

## To run
SpatialScanPoisson InputFile WindowInc WindowCount NumberofClusters NumberofMonteCarlo HighOrLowIndicator
 1. InputFile: a CSV without header that has 4 column: x, y, numberOfCases and intensity for each location
 2. WindowInc: the increment of scan windows
 3. WindowCount: the number of scan windows put at each location, and thus the maximum scan window is (WindowInc * WindowCount)
 4. NumberofClusters: number of clusters to detect
 5. NumberofMonteCarlo: number of Monte Carlo simulation
 6. HighOrLowIndicator:
	* 1: high value clusters only
	* -1: low value clusters only
	* 0: both
