# -*- Mode: Makefile; -*-

CCPLUS=g++
MPICC = mpic++
CFLAGS= -O3 -Wall -DTIMING
BINS=euler

all: $(BINS)

euler: *.cpp
	$(CCPLUS) $(CFLAGS) $^ -o $@

clean:
	 rm -f $(BINS) *.dat *.png