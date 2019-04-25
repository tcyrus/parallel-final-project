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
typedef enum { B, F, S, E, I, R, D, W } state;
const char stateNames[] = { 'B', 'F', 'S', 'E', 'I', 'R', 'D', 'W' }; // In case we actually want to print out the letters instead of numbers


/***************************************************************************/
/* Global Vars *************************************************************/
/***************************************************************************/

double g_time_in_secs = 0,
io_time_in_secs = 0;

unsigned long long g_start_cycles = 0,
g_end_cycles = 0,
io_start_cycles = 0,
io_end_cycles = 0;

const double g_processor_frequency = PROC_FREQ;
const double threshold = 0.75;

typedef struct {
	unsigned int time_in_state;
	state my_state;
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
	person->age = (unsigned int) ((rand() % 89) + 1);
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
			//Init based on population of board
			int chance = rand() % 100;
			Person *newPerson;
			newPerson = malloc(sizeof(Person));
			InitPerson(newPerson);

            if (i == 0 || j == 0 || i == height-1 || j == width-1) {
                newPerson->my_state = B;
            }
			else if (chance > POPULATION_RATE) {
				newPerson->my_state = S;
			}
			else {
				newPerson->my_state = F;
			}
			b->current[i][j] = *newPerson;

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

// Prints the board... Debug purposes 
void PrintBoard(Board* b) {
	for (size_t i = 0; i < b->height; i++) {
		for (size_t j = 0; j < b->width; j++) {
			printf(" %c ", stateNames[b->current[i][j].my_state]);
		}
		putchar('\n');
	}
}

/*
 Returns the number of people in a given state
 Takes a value from enum state as second argument
*/
unsigned int count_people(Board* b, state s) {
	unsigned int count = 0;
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
void infect_people(Board* b) {
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
	// TODO May change depending on if we keep border or not, probably should add the safe edge guarding regardless
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
*  TODO Should this function return a state or modify the state? probably return
*/
state next_state(Board* b, int x, int y) {
    Person* current_person;
    current_person = &b->current[x][y];
    state current_state = current_person->my_state;
	//Check current persons state to decide action
	if (current_state == S) { // If suceptible, count infected neighbors, decide if exposed
        int infectedNeighbors = get_infected_neighbors(b, x, y);
        int infectChance = (rand()%20) * infectedNeighbors;
        if (infectChance > 60) {
            return E;
        }
        else {
            return S;
        }
    }
    if (current_state == E) { // Decide if moves to Infected
        // TODO Need to decide on a incubation period for the infection to start
        if (current_person->time_in_state >= 5) {
            return I;
        }
        else{
        }
    }
    if (current_state == I) { // Decide if moves to W or D or R
        //TODO Time will decide if => W, other will decide if D or R
        if (current_person->time_in_state < 8) {
            int change = rand() % 100;
            if (change < 10) {
                int recovery = rand() % 100;
                if (recovery > 8) {
                    return R;
                }
                else {
                    return D;
                }
            }
            return I;
        }
        else{
            return W;
        }
    }
    if (current_state == W) { // Decide if moves to R or D

    }
    if (current_state == F) { // Decide if person will move into free cell

    }
    //If state is D, R, or B, will not change
    return current_state;
}
/*
 * Thread running function to perform the work for each tick on the cells
 */
void* rowTick(void* argp) {
	size_t num_rows = bc.height;
#if NUM_THREADS
	num_rows /= NUM_THREADS;
#endif
	size_t j = (*((size_t*)argp)) * num_rows;
	size_t end = j + num_rows - 1;
	for (; j < end; j++) {
		for (size_t i = 0; i < bc.width; i++) {
		    state current = bc.current[j][i].my_state;
		    state next = next_state(&bc, j, i);

		    if (current != next) {
		        bc.next[j][i].time_in_state = 0;
		        bc.next[j][i].my_state = next;
		    }
		    else {
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

/***************************************************************************/
/* Function: Main **********************************************************/
/***************************************************************************/

int main(int argc, char *argv[]) {
	g_start_cycles = (unsigned long long)GetTimeBase();

	// Using built in random for now, may change out later
	srand(time(NULL));

	InitBoard(&bc, GRID_SIZE+2, GRID_SIZE+2);

	unsigned int people = count_people(&bc, S);
	printf("Number of people: %d out of %d\n", people, GRID_SIZE*GRID_SIZE);

	infect_people(&bc);

	unsigned int infected = count_people(&bc, I);
	printf("Number of people infected: %d out of %d\n", infected, people);


	//Begin actual experiment
	for (int i=0; i<30; i++ ){
	    tick(&bc);
        infected = count_people(&bc, I) + count_people(&bc, W);
        printf("Number of people infected: %d out of %d\n", infected, people);
	}

    int infect_no_spread = count_people(&bc, W);
    printf("Number of people in W: %d out of %d\n", infect_no_spread, people);

    int recovered = count_people(&bc, R);
    printf("Number of people recovered: %d out of %d\n", recovered, people);

    int dead = count_people(&bc, D);
    printf("Number of people dead: %d out of %d\n", dead, people);


	g_end_cycles = (unsigned long long)GetTimeBase();
	g_time_in_secs = g_end_cycles - g_start_cycles;
	printf("Time = %f\n", g_time_in_secs);

	return 0;
}
