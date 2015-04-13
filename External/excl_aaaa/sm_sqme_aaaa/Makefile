####################################################
# M. Saimpert-interface for sm_aaaa coupling 
####################################################

.PHONY: clear test
.SUFFIXES: .cpp .o

CC=g++
CFLAGS=-g -Wall -DNUM_DOUBLE 
SOFLAGS=-fPIC
GSL=`gsl-config --cflags --libs`

all: sm_sqme_aaaa.a libsm_sqme_aaaa.so

objects=functions.o helicity_amplitudes.o sqme.o
f_0_objects=$(objects:%=src/%)

sm_sqme_aaaa.a:$(f_0_objects)
	ar -r $@ $(f_0_objects)

libsm_sqme_aaaa.so:$(f_0_objects)
	$(CC) $(SOFLAGS) -shared $(f_0_objects) -o $@ $(GSL)

.cpp.o:
	$(CC) $(CFLAGS) $(SOFLAGS) $(INS) -c $< -o $@ $(GSL)

clean:
	rm -rf src/*.o test/*.o *.a *.so testcall testcall2

