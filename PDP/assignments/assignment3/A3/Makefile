###############################################################################
# Makefile for assignment 3, Parallel and Distributed Computing.
###############################################################################

CC = mpicc
CFLAGS = -std=c99 -g -O3
LIBS = -lm

BIN = quicksort

all: $(BIN)

quicksort: quicksort.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
	
clean:
	$(RM) $(BIN)
