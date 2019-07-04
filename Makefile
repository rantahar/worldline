
CC=gcc 
CFLAGS= -march=native -Wall -Wextra -std=c99 -O3
LIB= -lm

DEPS=Makefile worldline.h mersenne.h

default: Thirring Thirring_exp

tests/test_worldline: mersenne_inline.o tests/test_worldline.c worldline.c
	$(CC) $(CFLAGS) -DTESTING tests/test_worldline.c worldline.c mersenne_inline.o -o tests/test_worldline -lcmocka -lm

.PHONY: tests
tests: tests/test_worldline
	tests/test_worldline

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

montecarlo: montecarlo.o worldline.o mersenne_inline.o $(DEPS)
	$(CC) $(CFLAGS) -o montecarlo montecarlo.o worldline.o mersenne_inline.o $(LIB)

LLR: LLR.o worldline.o mersenne_inline.o $(DEPS)
	$(CC) $(CFLAGS) -o LLR LLR.o worldline.o mersenne_inline.o $(LIB)

wanglandau: wanglandau.o worldline.o mersenne_inline.o  $(DEPS)
	$(CC) $(CFLAGS) -o wanglandau wanglandau.o worldline.o mersenne_inline.o $(LIB)

measure_sector: measure_sector.o worldline.o mersenne_inline.o  $(DEPS)
	$(CC) $(CFLAGS) -o measure_sector measure_sector.o worldline.o mersenne_inline.o $(LIB)

gauged: gauged.o worldline.o mersenne_inline.o $(DEPS)
	$(CC) $(CFLAGS) -o gauged gauged.o worldline.o mersenne_inline.o $(LIB)

clean:
	rm -f *.o 
	rm -f montecarlo LLR wanglandau measure_sector
