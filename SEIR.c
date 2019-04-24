/***************************************************************************/
/* Final Project	 ********************************************/
/* Timothy Cyrus and Liam Donohoe ******************************************/
/***************************************************************************/

/***************************************************************************/
/* Includes ****************************************************************/
/***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <stdbool.h>
#include <mpi.h>
#include <pthread.h>
#include <time.h>

//#define BGQ 1 // when running BG/Q, comment out when testing on mastiff

#ifdef BGQ
#include <hwi/include/bqc/A2_inlines.h>
#define PROC_FREQ 1600000000.0
#else
#define GetTimeBase MPI_Wtime
#define PROC_FREQ 1.0
#endif

/***************************************************************************/
/* Defines *****************************************************************/
/***************************************************************************/

#define GRID_SIZE 30
#define THREADS_PER_RANK 0
#define MAX_TICKS 256 // May look into making ticks represent hours, up to a number of days?
#define POPULATION_RATE 50 // Out of 100
#define INFECTION_RATE 10 // Out of 100
#define OUT_FILE "thresh.txt"

/*
	States include: SUBJECT TO CHANGE
	B - Border of automata
	F - Free Cell, no person, empty
	S - Suceptible to flu
	E - Exposed, incubation period, but not spread
	I - Infected, can spread flu
	R - Recovered, no longer spreads
	W - Infected without spreading
	D - Dead.
*/
typedef enum { B, F, S, E, I, R, D, W } state;
const char * stateNames[] = { "B", "F", "S", "E", "I", "R", "D", "W" }; // In case we actually want to print out the letters instead of numbers


/***************************************************************************/
/* Global Vars *************************************************************/
/***************************************************************************/

int world_size = -1,
world_rank = -1;

double g_time_in_secs = 0,
io_time_in_secs = 0;

unsigned long long g_start_cycles = 0,
g_end_cycles = 0,
io_start_cycles = 0,
io_end_cycles = 0;

const double g_processor_frequency = PROC_FREQ;

const double threshold = 0.75;

char * filename = OUT_FILE;

int* recv_above;
int* recv_below;
int* send_above;
int* send_below;

MPI_Status mpi_stat;
MPI_Request recieveUp,
recieveDown,
sendUp,
sendDown;

typedef struct {
	unsigned int time_in_state;
	state my_state;
	int age; //May enum? but need ranges of values so may not work

} Person;

typedef struct { //TODO Decide if using threads or ranks, if only threads, wont need width and height for grid
	size_t width;
	size_t height;
	Person** current;
	Person** previous;

} BoardChunk;

// May be a way to store the updated values before we finish going through all values
BoardChunk bc;

/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/


void Init_person(Person * person) {
	person->age = (rand() % 89) + 1;
	person->time_in_state = 0;
}

// Handles board initialization
// BoardChunk is clear on initialization
void Init_board(BoardChunk* b, size_t width, size_t height) {
	b->width = width;
	b->height = height;
	b->current = calloc(height, sizeof(Person*));
	b->previous = calloc(height, sizeof(Person*));

	for (size_t i = 0; i < height; i++) {

		b->current[i] = calloc(width, sizeof(Person));
		b->previous[i] = calloc(width, sizeof(Person));

		for (size_t j = 0; j < width; j++) {
			//Init based on population of board
			int chance = rand() % 100;
			Person *newPerson;
			newPerson = malloc(sizeof(Person));
			Init_person(newPerson);
			
			if (chance > POPULATION_RATE) {
				newPerson->my_state = S;
			}
			else {
				newPerson->my_state = F;
			}
			b->current[i][j] = *newPerson;

		}
	}
}

void Destroy_board(BoardChunk* b) {
	for (size_t i = 0; i < b->height; i++) {
		free(b->current[i]);
		free(b->previous[i]);
	}
	free(b->current);
	free(b->previous);
}

void Send_Recv_Init(size_t width) {
	// Each is one row of the current
	recv_above = (int*)calloc(width, sizeof(int));
	recv_below = (int*)calloc(width, sizeof(int));
	//send_above = (int*)calloc(width, sizeof(int));
	//send_below = (int*)calloc(width, sizeof(int));
}

void Send_Recv_Destroy() {
	// Each is one row of the current
	free(recv_above);
	free(recv_below);
	//free(send_above);
	//free(send_below);
}

// Prints the board... Debug purposes 
void print_board(BoardChunk* b) {
	for (size_t i = 0; i < b->height; i++) {
		for (size_t j = 0; j < b->width; j++) {
			printf(" %s ", stateNames[b->current[i][j].my_state]);
		}
		putchar('\n');
	}
}

/*
 Returns the number of people in a given state
 Takes a value from enum state as second argument
*/
int count_people(BoardChunk* b, state s) {
	int count = 0;
	for (size_t i = 0; i < b->height; i++) {
		for (size_t j = 0; j < b->width; j++) {
			if (b->current[i][j].my_state == s)
				count++;
		}
	}
	return count;
}

/*
* Uses the define of INFECTION_RATE to randomly infect a percent of the population
* This is only to be used as an initialization function before beginning actual testing
*/
void infect_people(BoardChunk* b) {
	//We will infect a certain percent of the population as a starting point
	for (size_t i = 0; i < b->height; i++) {
		for (size_t j = 0; j < b->width; j++) {
			if (b->current[i][j].my_state != S) {
				int infected = rand() % 100;
				if (infected < INFECTION_RATE) {
					b->current[i][j].my_state = I;
				}
			}	
		}
	}
}

/*
* BoardChunk function to obtain the number of infected neighbors any given cell has
* Parameters: x and y coordinates of target cell
* Output: number of infected cells adjacent to the target cell
* We will make the border of the world be a border, so this function will 
*	never recieve a cell on the edge of the grid
*/
size_t get_infected_neighbors(BoardChunk* b, int x, int y) {
	//We will not be using wrap around for this version
	size_t count = 0;
	int left = x - 1,
		right = x + 1,
		up = y + 1,
		down = y - 1;
	

	//We only care about infected cells in state I,
	//	infected cells in state W are non-infectious
	if (b->current[left][y].my_state == I)
		count++;
	if (b->current[left][up].my_state == I)
		count++;
	if (b->current[left][down].my_state == I)
		count++;
	if (b->current[x][up].my_state == I)
		count++;
	if (b->current[x][down].my_state == I)
		count++;
	if (b->current[right][y].my_state == I)
		count++;
	if (b->current[right][up].my_state == I)
		count++;
	if (b->current[right][down].my_state == I)
		count++;

	return count;

}

/*
* This function will compute the next state for any given cell at the current time
* Parameters: x and y coordinates of target cell, and pointer to board itself
* Output: The value of the state that will succeed this state
* 
*/
void next_state(BoardChunk * b, int x, int y) {
	Person * current_person = &b->current[x][y];
	//Check current persons state to decide action
	if (current_person->my_state == S) { // If suceptible, count infected neighbors, decide if exposed
		int infectedNeighbors = get_infected_neighbors(b, x, y);
		printf("Found %d infected neigbors\n", infectedNeighbors);
		if (infectedNeighbors > 0) {
			current_person->my_state = E;
		}
	}


}

void* rowTick(void* argp) {
	size_t num_rows = bc.height;
#if THREADS_PER_RANK
	num_rows /= THREADS_PER_RANK;
#endif
	size_t j = (*((size_t*)argp)) * num_rows;
	int end = j + num_rows - 1;

	for (; j < end; j++) {
		for (size_t i = 0; i < bc.width; i++) {

		}
	}

	return NULL;
}
/*
void tick(BoardChunk* b) {
	int rankUp = (world_rank == 0) ? (world_size - 1) : (world_rank - 1),
		rankDown = (world_rank == (world_size - 1)) ? 0 : (world_rank + 1);

	//memcpy(send_above, b->current[0], sizeof(b->current));
	send_above = b->current[0];
	//memcpy(send_below, b->current[b->height-1], sizeof(b->current));
	send_below = b->current[b->height - 1];

	for (size_t i = 0; i < b->height; i++) {
		memcpy(b->previous[i], b->current[i], sizeof(b->current));
	}

	// Send at the beginning of the tick
	// Send top of chunk to rank+1
	MPI_Isend(send_above, b->width, MPI_INT, rankUp, 0, MPI_COMM_WORLD, &sendUp);

	// Send bottom of chunk to rank-1
	MPI_Isend(send_below, b->width, MPI_INT, rankDown, 0, MPI_COMM_WORLD, &sendDown);

	// Recieve from rank-1 as bottom
	MPI_Irecv(recv_above, b->width, MPI_INT, rankUp, 0, MPI_COMM_WORLD, &recieveUp);

	// Recieve from rank+1 as top
	MPI_Irecv(recv_below, b->width, MPI_INT, rankDown, 0, MPI_COMM_WORLD, &recieveDown);

#if THREADS_PER_RANK
	pthread_t* tid = calloc(THREADS_PER_RANK, sizeof(pthread_t));
#endif
	size_t* tmpI = calloc(THREADS_PER_RANK, sizeof(size_t));

	MPI_Wait(&recieveUp, &mpi_stat);
	MPI_Wait(&recieveDown, &mpi_stat);

	pthread_t parent = pthread_self();

	tmpI[0] = 0;

#if THREADS_PER_RANK
	// Also need to keep track of current board without changes
	for (size_t i = 1; i < THREADS_PER_RANK; i++) {
		tmpI[i] = i;
		pthread_create(&(tid[i]), NULL, rowTick, (void*)&tmpI[i]);
	}
#endif

	MPI_Wait(&sendUp, &mpi_stat);
	MPI_Wait(&sendDown, &mpi_stat);

	if (pthread_self() == parent) {
		rowTick((void*)&tmpI[0]);
	}

#if THREADS_PER_RANK
	for (size_t i = 1; i < THREADS_PER_RANK; i++) {
		pthread_join(tid[i], NULL);
	}

	free(tid);
#endif
	free(tmpI);

	int** tmp = b->current;
	b->current = b->previous;
	b->previous = tmp;
}*/

/***************************************************************************/
/* Function: Main **********************************************************/
/***************************************************************************/

int main(int argc, char *argv[]) {
	// Example MPI startup and using CLCG4 RNG
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);


	if (world_rank == 0) {
		g_start_cycles = (unsigned long long)GetTimeBase();
	}
	//Using built in random for now, may change out later
	srand(time(NULL));

	size_t chunk = GRID_SIZE / world_size;
	Init_board(&bc, GRID_SIZE, chunk);
	Send_Recv_Init(bc.width);

	print_board(&bc);

	int people = count_people(&bc, S);
	printf("Number of people: %d out of %d\n", people, GRID_SIZE*GRID_SIZE);

	infect_people(&bc);

	int infected = count_people(&bc, I);
	printf("Number of people infected: %d out of %d\n", infected, people);

	print_board(&bc);


	if (world_rank == 0) {
		g_end_cycles = (unsigned long long)GetTimeBase();
		g_time_in_secs = g_end_cycles - g_start_cycles;
		printf("Time = %f\n", g_time_in_secs);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	// END -Perform a barrier and then leave MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}
