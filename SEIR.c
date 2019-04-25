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
#define NUM_THREADS 4
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
	Person** next;
} Board;

// May be a way to store the updated values before we finish going through all values
Board bc;

/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

void InitPerson(Person* person) {
	person->age = (unsigned int) ((random() % 89) + 1);
	person->time_in_state = 0;
}

// Handles board initialization
// BoardChunk is clear on initialization
void InitBoard(Board* b, size_t width, size_t height) {
	b->width = width;
	b->height = height;
	b->current = calloc(height, sizeof(Person*));
	b->next = calloc(height, sizeof(Person*));

	for (size_t i = 0; i < height; i++) {
		b->current[i] = calloc(width, sizeof(Person));
		b->next[i] = calloc(width, sizeof(Person));

		for (size_t j = 0; j < width; j++) {
			// Init based on population of board
			int chance = random() % 100;
			InitPerson(&(b->current[i][j]));

            if (i == 0 || j == 0 || i == height-1 || j == width-1) {
                b->current[i][j].state = BOARDER_CELL;
            } else if (chance > POPULATION_RATE) {
				b->current[i][j].state = SUCEPTIBLE_CELL;
			} else {
				b->current[i][j].state = FREE_CELL;
			}
		}
	}
}

void DestroyBoard(Board* b) {
	for (size_t i = 0; i < b->height; i++) {
		free(b->current[i]);
		free(b->next[i]);
	}
	free(b->current);
	free(b->next);
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
			if (b->current[i][j].state != SUCEPTIBLE_CELL) {
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
int get_infected_neighbors(Board* b, int x, int y) {
	// We will not be using wrap around for this version
	size_t count = 0;
	int left = x - 1,
		right = x + 1,
		up = y + 1,
		down = y - 1;

	// We only care about infected cells in state I,
	// infected cells in state W are non-infectious
    count += (b->current[x][up].state == INFECTED_CELL);
    count += (b->current[x][down].state == INFECTED_CELL);
	if (left >= 0) {
        count += (b->current[left][y].state == INFECTED_CELL);
        count += (y == b->height - 1) && (b->current[left][up].state == INFECTED_CELL);
        count += (down >= 0) && (b->current[left][down].state == INFECTED_CELL);
    }
    if (right < b->width) {
        count += (b->current[right][y].state == INFECTED_CELL);
        count += (y == b->height - 1) && (b->current[right][up].state == INFECTED_CELL);
        count += (down >= 0) && (b->current[right][down].state == INFECTED_CELL);
    }

	return count;

}

/*
* This function will compute the next state for any given cell at the current time
* Parameters: x and y coordinates of target cell, and pointer to board itself
* Output: The value of the state that will succeed this state
*  TODO Should this function return a state or modify the state? probably return
*/
cell_state next_state(Board* b, int x, int y) {
    Person* current_person = &(b->current[x][y]);

	// Check current persons state to decide action
	switch (current_person->state) {
	    case SUCEPTIBLE_CELL: ;
            // If suceptible, count infected neighbors, decide if exposed
            int infectedNeighbors = get_infected_neighbors(b, x, y);
            int infectChance = (random() % 20) * infectedNeighbors;
            return (infectChance > 60) ? EXPOSED_CELL : SUCEPTIBLE_CELL;
        case EXPOSED_CELL:
            // Decide if moves to Infected
            // TODO Need to decide on a incubation period for the infection to start
            if (current_person->time_in_state >= 5) {
                return INFECTED_CELL;
            } else {}
            break;
	    case INFECTED_CELL:
            // Decide if moves to W or D or R
            //TODO Time will decide if => W, other will decide if D or R
            if (current_person->time_in_state < 8) {
                int change = random() % 100;
                if (change < 10) {
                    int recovery = random() % 100;
                    return (recovery > 8) ? RECOVERED_CELL : DEAD_CELL;
                }

                return INFECTED_CELL;
            }

            return WITHOUT_CELL;
        case WITHOUT_CELL:
            // Decide if moves to R or D
            break;
	    case FREE_CELL:
	        // Decide if person will move into free cell
	        break;
	    default:
	        break;
    }

    // If state is D, R, or B, will not change
    return current_person->state;
}

/*
 * Thread running function to perform the work for each tick on the cells
 */
void* rowTick(void* argp) {
    unsigned int num_rows = bc.height;
#if NUM_THREADS
	num_rows /= NUM_THREADS;
#endif
	unsigned int j = (*((size_t*)argp)) * num_rows;
	unsigned int end = j + num_rows - 1;
	for (; j < end; j++) {
		for (size_t i = 0; i < bc.width; i++) {
		    cell_state current = bc.current[j][i].state;
		    cell_state next = next_state(&bc, j, i);

		    if (current != next) {
		        bc.next[j][i].time_in_state = 0;
		        bc.next[j][i].state = next;
		    } else {
                bc.next[j][i].time_in_state++;
            }
		}
	}

	return NULL;
}

/*
 * Master function to handle all actions for each time unit
 */
void tick(Board* b) {
	// Copy Board
	for (size_t i = 0; i < b->height; i++) {
		memcpy(b->next[i], b->current[i], sizeof(b->current));
	}

#if NUM_THREADS
	pthread_t* tid = calloc(NUM_THREADS, sizeof(pthread_t));
#endif
	size_t* tmpI = calloc(NUM_THREADS, sizeof(size_t));

	pthread_t parent = pthread_self();

	tmpI[0] = 0;

#if NUM_THREADS
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
	b->current = b->next;
	b->next = tmp;
}


int main() {
	g_start_cycles = (unsigned long long)GetTimeBase();

	// Using built in random for now, may change out later
	srandom(time(NULL));

	InitBoard(&bc, GRID_SIZE+2, GRID_SIZE+2);

    infect_people(&bc);

    unsigned int people, infected;

    people = count_people(&bc, SUCEPTIBLE_CELL);
	printf("Number of people suceptible: %d out of %d\n", people, GRID_SIZE*GRID_SIZE);

	infected = count_people(&bc, INFECTED_CELL);
	printf("Number of people infected: %d out of %d\n", infected, people);

#ifdef DEBUG
    PrintBoard(&bc);
#endif

	// Begin actual experiment
	for (size_t i = 0; i < 30; i++) {
	    tick(&bc);
        infected = count_people(&bc, INFECTED_CELL) + count_people(&bc, WITHOUT_CELL);
        printf("Number of people infected: %d out of %d\n", infected, people);
	}

    int infect_no_spread = count_people(&bc, WITHOUT_CELL);
    printf("Number of people in W: %d out of %d\n", infect_no_spread, people);

    int recovered = count_people(&bc, RECOVERED_CELL);
    printf("Number of people recovered: %d out of %d\n", recovered, people);

    int dead = count_people(&bc, DEAD_CELL);
    printf("Number of people dead: %d out of %d\n", dead, people);

#ifdef DEBUG
    PrintBoard(&bc);
#endif

    DestroyBoard(&bc);

    g_end_cycles = (unsigned long long)GetTimeBase();
	g_time_in_secs = g_end_cycles - g_start_cycles;
	printf("Time = %f\n", g_time_in_secs);

	return 0;
}
