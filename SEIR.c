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
#include <pthread.h>
#include <time.h>

#ifdef __bgq__
#include <hwi/include/bqc/A2_inlines.h>
#define PROC_FREQ 1600000000.0
#else
double GetTimeBase() {
	return (double)clock() / CLOCKS_PER_SEC;
}
#define PROC_FREQ 1.0
#endif

/***************************************************************************/
/* Defines *****************************************************************/
/***************************************************************************/

#define GRID_SIZE 30
#define NUM_THREADS 2
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
typedef enum {
    BOARDER_CELL,
    FREE_CELL,
    SUCEPTIBLE_CELL,
    EXPOSED_CELL,
    INFECTED_CELL,
    RECOVERED_CELL,
    DEAD_CELL,
    WITHOUT_CELL
} cell_state;

#ifdef DEBUG
// In case we actually want to print out the letters instead of numbers
const char stateNames[] = { 'B', 'F', 'S', 'E', 'I', 'R', 'D', 'W' };
#endif


/***************************************************************************/
/* Global Vars *************************************************************/
/***************************************************************************/

double g_time_in_secs = 0;
double io_time_in_secs = 0;

unsigned long long g_start_cycles = 0;
unsigned long long g_end_cycles = 0;
unsigned long long io_start_cycles = 0;
unsigned long long io_end_cycles = 0;

const double g_processor_frequency = PROC_FREQ;
const double threshold = 0.75;

typedef struct {
	unsigned int time_in_state;
    cell_state state;
	unsigned int age;
} Person;

typedef struct {
	size_t width;
	size_t height;
	Person** current;
	Person** previous;
} Board;

// May be a way to store the updated values before we finish going through all values
Board bc;

/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

void InitPerson(Person* person) {
	person->age = (random() % 89) + 1;
	person->time_in_state = 0;
}

// Handles board initialization
// BoardChunk is clear on initialization
void InitBoard(Board* b, size_t width, size_t height) {
	b->width = width;
	b->height = height;
	b->current = calloc(height, sizeof(Person*));
	b->previous = calloc(height, sizeof(Person*));

	for (size_t i = 0; i < height; i++) {
		b->current[i] = calloc(width, sizeof(Person));
		b->previous[i] = calloc(width, sizeof(Person));

		for (size_t j = 0; j < width; j++) {
			//Init based on population of board
			int chance = random() % 100;
			InitPerson(&(b->current[i][j]));

            b->current[i][j].state = (chance > POPULATION_RATE) ? SUCEPTIBLE_CELL : FREE_CELL;
		}
	}
}

void DestroyBoard(Board* b) {
	for (size_t i = 0; i < b->height; i++) {
		free(b->current[i]);
		free(b->previous[i]);
	}
	free(b->current);
	free(b->previous);
}

#ifdef DEBUG
// Prints the board... Debug purposes 
void PrintBoard(Board* b) {
	for (size_t i = 0; i < b->height; i++) {
		for (size_t j = 0; j < b->width; j++) {
			printf(" %c ", stateNames[b->current[i][j].state]);
		}
		putchar('\n');
	}
}
#endif

/*
 Returns the number of people in a given state
 Takes a value from enum state as second argument
*/
unsigned int count_people(Board* b, cell_state cell) {
	unsigned int count = 0;
	for (size_t i = 0; i < b->height; i++) {
		for (size_t j = 0; j < b->width; j++) {
			if (b->current[i][j].state == cell)
				count++;
		}
	}
	return count;
}

/*
* Uses the define of INFECTION_RATE to randomly infect a percent of the population
* This is only to be used as an initialization function before beginning actual testing
*/
void infect_people(Board* b) {
	//We will infect a certain percent of the population as a starting point
	for (size_t i = 0; i < b->height; i++) {
		for (size_t j = 0; j < b->width; j++) {
			if (b->previous[i][j].state != SUCEPTIBLE_CELL) {
				int infected = random() % 100;
				if (infected < INFECTION_RATE) {
					b->current[i][j].state = INFECTED_CELL;
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
* never recieve a cell on the edge of the grid
*/
size_t get_infected_neighbors(Board* b, int x, int y) {
	// We will not be using wrap around for this version
	size_t count = 0;
	int left = x - 1,
		right = x + 1,
		up = y + 1,
		down = y - 1;

	// We only care about infected cells in state I,
	// infected cells in state W are non-infectious
    count += (b->previous[x][up].state == INFECTED_CELL);
    count += (b->previous[x][down].state == INFECTED_CELL);
	if (left >= 0) {
        count += (b->previous[left][y].state == INFECTED_CELL);
        count += (up < b->height) && (b->previous[left][up].state == INFECTED_CELL);
        count += (down >= 0) && (b->previous[left][down].state == INFECTED_CELL);
    }
    if (right < b->width) {
        count += (b->previous[right][y].state == INFECTED_CELL);
        count += (up < b->height) && (b->previous[right][up].state == INFECTED_CELL);
        count += (down >= 0) && (b->previous[right][down].state == INFECTED_CELL);
    }

	return count;

}

/*
* This function will compute the next state for any given cell at the current time
* Parameters: x and y coordinates of target cell, and pointer to board itself
* Output: The value of the state that will succeed this state
*/
void next_state(Board* b, int x, int y) {
	Person* p = &b->current[x][y];
	// Check current persons state to decide action
	if (p->state == SUCEPTIBLE_CELL) {
	    // If suceptible, count infected neighbors, decide if exposed
		int infectedNeighbors = get_infected_neighbors(b, x, y);
#ifdef DEBUG
		printf("Found %d infected neigbors\n", infectedNeighbors);
#endif
		if (infectedNeighbors > 0) {
			p->state = EXPOSED_CELL;
		}
	}
}

void* rowTick(void* argp) {
    unsigned int num_rows = bc.height;
#if NUM_THREADS
	num_rows /= NUM_THREADS;
#endif
    unsigned int j = (*((size_t*)argp)) * num_rows;
	unsigned int end = j + num_rows - 1;

	for (; j < end; j++) {
		for (size_t i = 0; i < bc.width; i++) {
            next_state(&bc, i, j);
		}
	}

	return NULL;
}

void tick(Board* b) {
	// Copy Board
	for (size_t i = 0; i < b->height; i++) {
		memcpy(b->previous[i], b->current[i], sizeof(b->current));
	}

#if NUM_THREADS
	pthread_t* tid = calloc(NUM_THREADS, sizeof(pthread_t));
#endif
	size_t* tmpI = calloc(NUM_THREADS, sizeof(size_t));

	pthread_t parent = pthread_self();

	tmpI[0] = 0;

#if NUM_THREADS
	// Also need to keep track of current board without changes
	for (size_t i = 1; i < NUM_THREADS; i++) {
		tmpI[i] = i;
		pthread_create(&(tid[i]), NULL, rowTick, (void*)&tmpI[i]);
	}
#endif

	if (pthread_self() == parent) {
		rowTick((void*)&tmpI[0]);
	}

#if NUM_THREADS
	for (size_t i = 1; i < NUM_THREADS; i++) {
		pthread_join(tid[i], NULL);
	}

	free(tid);
#endif
	free(tmpI);

	Person** tmp = b->current;
	b->current = b->previous;
	b->previous = tmp;
}


int main() {
	g_start_cycles = (unsigned long long)GetTimeBase();

	// Using built in random for now, may change out later
	srandom(time(NULL));

	InitBoard(&bc, GRID_SIZE, GRID_SIZE);

#ifdef DEBUG
	PrintBoard(&bc);
#endif

	unsigned int people = count_people(&bc, SUCEPTIBLE_CELL);
	printf("Number of people suceptible: %d out of %d\n", people, GRID_SIZE*GRID_SIZE);

	infect_people(&bc);

	unsigned int infected = count_people(&bc, INFECTED_CELL);
	printf("Number of people infected: %d out of %d\n", infected, people);

#ifdef DEBUG
    PrintBoard(&bc);
#endif

    for (size_t i = 0; i < MAX_TICKS; i++) {
        // Perform tick
        tick(&bc);
    }

#ifdef DEBUG
    PrintBoard(&bc);
#endif

    DestroyBoard(&bc);

    g_end_cycles = (unsigned long long)GetTimeBase();
	g_time_in_secs = g_end_cycles - g_start_cycles;
	printf("Time = %f\n", g_time_in_secs);

	return 0;
}
