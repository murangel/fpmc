####################################################
# Interface to anomalous aaww coupling from comphep
# O. Kepka 13/12/2008 kepkao@fzu.cz
# Original comphep rutines are in directory comphep
# Model generated files + call_util.c in src
####################################################

.PHONY: clear test
.SUFFIXES: .cpp .o

INS=-Icomphep/include

#CHEP=/home/olda/work/mc/comphep/install/comphep-4.4p3-gcc411
#INS:=-I$(CHEP)/include
#CC=    $(shell more $(CHEP)/CC)
#FRT=   $(shell more $(CHEP)/FRT)
#CFLAGS=$(shell more $(CHEP)/CFLAGS)
#CLIBS= $(shell more $(CHEP)/CLIBS) 

CC=g++
CFLAGS=-g -Wall -DNUM_DOUBLE -DCOMPHEP

all: sqme_aaww.a
test:testcall testcall2

objects=d1.o f1.o f2.o f3.o f4.o f5.o f6.o  service.o sqme.o call_util.o
f_0_objects=$(objects:%=src/%)

sqme_aaww.a:$(f_0_objects)
	ar -r $@ $(f_0_objects)

.cpp.o:
	$(CC) $(CFLAGS) $(INS) -c $< -o $@

testcall:test/testcall.o sqme_aaww.a
	$(CC) $(CFLAGS) -o $@ -lm test/testcall.o sqme_aaww.a 

testcall2:test/testcall2.o sqme_aaww.a
	$(CC) $(CFLAGS) -o $@ -lm test/testcall2.o sqme_aaww.a 
clean:
	rm -rf src/*.o test/*.o *.a testcall testcall2

