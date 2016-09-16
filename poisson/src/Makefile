GCC= g++
OMPFLAGS = -g -fopenmp

all: SpatialScanPoisson

io.o: io.c
	$(GCC) -o $@ -c $<
scan.o: scan.c
	$(GCC) $(OMPFLAGS) -o $@ -c $<
mc.o: mc.c
	$(GCC) $(OMPFLAGS) -o $@ -c $< -std=c++11
spatialScanPoisson.o: spatialScanPoisson.c
	$(GCC) -o $@ -c $<
SpatialScanPoisson: spatialScanPoisson.o io.o scan.o mc.o
	$(GCC) $(OMPFLAGS) -o ../$@ $+
clean:
	rm -f *.o ../SpatialScanPoisson

