Running SEIR.c

Our program can be run and compiled with or without MPI. Compiling is normal, all additional files are included in SEIR.c. A basic compilation line with MPI is included below for your convinience.

mpicc -I. -pthread -g -Wall -O5 SEIR.c -o SEIR.exe

All of the program arguments can be found at the top of the SEIR.c file, and include, 

GRID_SIZE 64 : The size of the Cellular Automata to be run
NUM_THREADS 4 : Number of threads to use for computing states
MAX_TICKS 1024 : Maximum running time for the Cellular Automata, will terminate if reaches this point
POPULATION_RATE 50 (Out of 100) : Percentage of the grid to populate with living person objects
INFECTION_RATE 1 (Out of 100) : Initial ratio to infect the grid with
RECOVERY_RATE 3 : Rate at which a person will recover at, based on number of days and this value
INFECTIVE_TIME 8 : Amount of time the person remains able to infect other people

When running SEIR.exe, all results will be printed whether using multiple threads or multiple ranks. 
