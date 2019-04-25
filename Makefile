CC      = gcc
MPICC   = mpicc
CFLAGS  = -g -Wall -O5
LDFLAGS = -I.

all: SEIR.exe

SEIR.exe: SEIR.c
	$(MPICC) $(LDFLAGS) -pthread $(CFLAGS) SEIR.c -o SEIR.exe

.PHONY: clean
clean:
	$(RM) SEIR.exe
