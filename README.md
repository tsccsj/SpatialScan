# Spatial Scan Statistics
This is an implementation of Martin Kulldorff's Spatial Scan Statistics(1997) for finding spatial clusters (local excesses of events) in a spatial point process.   
It is implemented using C with OpenMP for parallelization. 
Two models are implemented: a *Poisson Model* and a *Bernoulli Model*. 

## Poisson Model
A Poisson model deals with the number of events occurring in a time interval and a spatial region.   
### To run
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


## Bernoulli Model
A Bernoulli model handles events that are in either one of two states (i.e., belonging to either of two categories), which is often used to compare the spatial distributions of two types of events, such as a case control study.   
## To run
SpatialScanBernoulli InputFile WindowInc WindowCount NumberofClusters NumberofMonteCarlo HighOrLowIndicator
 1. InputFile: a CSV without header that has 4 column: x, y, numberOfCases and numberOfControls for each location
 2. WindowInc: the increment of scan windows
 3. WindowCount: the number of scan windows put at each location, and thus the maximum scan window is (WindowInc * WindowCount)
 4. NumberofClusters: number of clusters to detect
 5. NumberofMonteCarlo: number of Monte Carlo simulation
 6. HighOrLowIndicator:
	* 1: high value clusters only
	* -1: low value clusters only
	* 0: both
 

## Reference
Kulldorff, M., 1997. A spatial scan statistic. Communications in Statistics-Theory and methods, 26(6), pp.1481-1496.
