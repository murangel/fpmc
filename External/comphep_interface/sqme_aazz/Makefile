####################################################
# Interface to anomalous aazz coupling from comphep
# O. Kepka 21/04/2009 kepkao@fzu.cz
# Original comphep rutines are in directory comphep
# Model generated files + call_util.c in src
####################################################

.PHONY: clear test
.SUFFIXES: .cpp .o

INS=-Icomphep/include

CC=g++
CFLAGS=-g -Wall -DNUM_DOUBLE -DCOMPHEP

all: sqme_aazz.a
test:testcall testcall2

objects=d1.o f1.o service.o sqme.o call_util.o
f_0_objects=$(objects:%=src/%)

sqme_aazz.a:$(f_0_objects)
	ar -r $@ $(f_0_objects)

.cpp.o:
	$(CC) $(CFLAGS) $(INS) -c $< -o $@

testcall:test/testcall.o sqme_aazz.a
	$(CC) $(CFLAGS) -o $@ -lm test/testcall.o sqme_aazz.a 

testcall2:test/testcall2.o sqme_aazz.a
	$(CC) $(CFLAGS) -o $@ -lm test/testcall2.o sqme_aazz.a 

clean:
	rm -rf src/*.o test/*.o *.a testcall testcall2

