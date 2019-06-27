
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

worldline: worldline.o mersenne_inline.o $(DEPS)
	$(CC) $(CFLAGS) -o worldline worldline.o mersenne_inline.o $(LIB)

worldline_llr.o: worldline.c $(DEPS)
	$(CC) $(CFLAGS) -DLLR -c -o $@ $< $(CFLAGS)

LLR: LLR.o mersenne_inline.o worldline_llr.o  $(DEPS)
	$(CC) $(CFLAGS) -o LLR LLR.o worldline_llr.o mersenne_inline.o $(LIB)

wanglandau: wanglandau.o worldline.o mersenne_inline.o  $(DEPS)
	$(CC) $(CFLAGS) -o wanglandau wanglandau.o worldline.o mersenne_inline.o $(LIB)

worldline_sector.o: worldline.c $(DEPS)
	$(CC) $(CFLAGS) -DMEASURE_SECTOR -c -o $@ $< $(CFLAGS)

worldline_sector: worldline_sector.o mersenne_inline.o  $(DEPS)
	$(CC) $(CFLAGS) -DMEASURE_SECTOR -o worldline_sector worldline_sector.o mersenne_inline.o $(LIB)


clean:
	rm -f *.o  
