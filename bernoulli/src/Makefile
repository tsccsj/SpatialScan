GCC= g++
OMPFLAGS = -g -fopenmp

all: SpatialScanBernoulli

io.o: io.c
	$(GCC) -o $@ -c $<
scan.o: scan.c
	$(GCC) $(OMPFLAGS) -o $@ -c $<
mc.o: mc.c
	$(GCC) $(OMPFLAGS) -o $@ -c $< -std=c++11
spatialScanBernoulli.o: spatialScanBernoulli.c
	$(GCC) -o $@ -c $<
SpatialScanBernoulli: spatialScanBernoulli.o io.o scan.o mc.o
	$(GCC) $(OMPFLAGS) -o ../$@ $+
clean:
	rm -f *.o ../SpatialScanBernoulli

