###############################################################################
# Makefile for project, Parallel and Distributed Computing.
###############################################################################

CC = mpicc
CFLAGS = -std=c99 -g -O3 -Wall -S
LIBS = -lm

BIN = project

all: $(BIN)

project: project.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
	
clean:
	$(RM) $(BIN)
