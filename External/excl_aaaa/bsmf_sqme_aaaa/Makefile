####################################################
# M. Saimpert-interface for bsmf_aaaa coupling 
####################################################

.PHONY: clear test
.SUFFIXES: .cpp .o

CC=g++
CFLAGS=-g -Wall -DNUM_DOUBLE 
SOFLAGS=-fPIC
GSL=`gsl-config --cflags --libs`

all: bsmf_sqme_aaaa.a libbsmf_sqme_aaaa.so

objects=functions.o helicity_amplitudes.o sqme.o
f_0_objects=$(objects:%=src/%)

bsmf_sqme_aaaa.a:$(f_0_objects)
	ar -r $@ $(f_0_objects)

libbsmf_sqme_aaaa.so:$(f_0_objects)
	$(CC) $(SOFLAGS) -shared $(f_0_objects) -o $@ $(GSL)

.cpp.o:
	$(CC) $(CFLAGS) $(SOFLAGS) $(INS) -c $< -o $@ $(GSL)

clean:
	rm -rf src/*.o test/*.o *.a *.so testcall testcall2

