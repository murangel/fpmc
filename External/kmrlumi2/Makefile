# Makefile for MSTW Fortran and C++ code.
# Comments to Graeme Watt <watt(at)hep.ucl.ac.uk>.

FC=g77 -ggdb
LDLIBS=-lg2c 
#FC=gfortran
#LDLIBS=-lgfortran
CC=g++
FFLAGS=-Wall -W
CCFLAGS=-Wall -W
EXEC=GPDF

all: $(EXEC)

#####################################
# Fortran version:

GPDF : GPDF.o mstwpdf.o alphaS.o 
	$(FC) $^ -o $@

GPDF.o : GPDF.f
	$(FC) -c $(FFLAGS) $< -o $@

mstwpdf.o : mstwpdf.f
	$(FC) -c $(FFLAGS) $< -o $@

alphaS.o : alphaS.f
	$(FC) -c $(FFLAGS) $< -o $@

clean :
	-rm -f *.o $(EXEC)
