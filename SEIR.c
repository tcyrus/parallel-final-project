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
#include <time.h>

#include <mpi.h>
#include <pthread.h>

#include "cell_state.h"
#include "person.h"
#include "board.h"

#ifdef __bgq__
#include <hwi/include/bqc/A2_inlines.h>
#define PROC_FREQ 1600000000.0
#else
#define GetTimeBase MPI_Wtime
/*double GetTimeBase() { return (double)clock() / CLOCKS_PER_SEC; }*/
#define PROC_FREQ 1.0
#endif

/***************************************************************************/
/* Defines *****************************************************************/
/***************************************************************************/

#define GRID_SIZE 30
#define NUM_THREADS 8
// May look into making ticks represent hours, up to a number of days?
#define MAX_TICKS 256
#define POPULATION_RATE 50 // Out of 100
#define INFECTION_RATE 10 // Out of 100
#define OUT_FILE "thresh.txt"

Person* recv_above;
Person* recv_below;


/***************************************************************************/
/* Global Vars *************************************************************/
/***************************************************************************/

int world_size = -1;
int world_rank = -1;

double g_time_in_secs = 0;

unsigned long long g_start_cycles = 0;
unsigned long long g_end_cycles = 0;

Person* send_above;
Person* send_below;
Person* recv_above;
Person* recv_below;

MPI_Status mpi_stat;
MPI_Request recieveUp,
            recieveDown,
            sendUp,
            sendDown;

// May be a way to store the updated values before we finish going through all values
Board bc;

/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

void Send_Recv_Init(size_t width) {
    // Each is one row of the grid
    recv_above = (Person*)calloc(width, sizeof(int));
    recv_below = (Person*)calloc(width, sizeof(int));
}

void Send_Recv_Destroy() {
    // Each is one row of the grid
    free(recv_above);
    free(recv_below);
}

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

            /*if (i == 0 || j == 0 || i == height-1 || j == width-1) {
                b->current[i][j].state = BOARDER_CELL;
            } else */if (chance > POPULATION_RATE) {
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
			printf(" %c ", cell_state_names[b->current[i][j].state]);
		}
		putchar('\n');
	}
    putchar('\n');
}
#endif

/*
 Returns the number of people in a given state
 Takes a value from enum state as second argument
*/
unsigned long long count_people(Board* b, cell_state cell) {
	unsigned long long count = 0;
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
	// We will infect a certain percent of the population as a starting point
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
size_t get_infected_neighbors(Board* b, unsigned int x, unsigned int y) {
	// We will not be using wrap around for this version
	size_t count = 0;
    int left = (x == 0) ? (b->width - 1) : (x - 1);
    int right = (x == b->width - 1) ? 0 : (x + 1);

	// We only care about infected cells in state I,
	// infected cells in state W are non-infectious
    if (y == 0) {
        count += (recv_below[x].state == INFECTED_CELL); // down 1
        count += (recv_below[left].state == INFECTED_CELL); // Diagonal left down -> wrap to right edge
        count += (recv_below[right].state == INFECTED_CELL); // Diagonal right down
    } else {
        count += (b->current[y-1][right].state == INFECTED_CELL); // down and right
        count += (b->current[y-1][x].state == INFECTED_CELL); // just down
        count += (b->current[y-1][left].state == INFECTED_CELL); // down and left
    }


    if (y == b->height - 1) {
        count += (recv_above[x].state == INFECTED_CELL); // up 1
        count += (recv_above[left].state == INFECTED_CELL); // Diagonal left up -> wrap to right edge
        count += (recv_above[right].state == INFECTED_CELL); // Diagonal right up
    } else {
        count += (b->current[y+1][x].state == INFECTED_CELL); // Just up
        count += (b->current[y+1][left].state == INFECTED_CELL); // UP and left
        count += (b->current[y+1][right].state == INFECTED_CELL); // UP and right
    }

    count += (b->current[y][right].state == INFECTED_CELL); // just right
    count += (b->current[y][left].state == INFECTED_CELL); // just left

	return count;

}

/*
* This function will compute the next state for any given cell at the current time
* Parameters: x and y coordinates of target cell, and pointer to board itself
* Output: The value of the state that will succeed this state
*  TODO Should this function return a state or modify the state? probably return
*/
cell_state next_state(Board* b, unsigned int x, unsigned int y) {
    Person* current_person = &(b->current[x][y]);

	// Check current persons state to decide action
	switch (current_person->state) {
	    case SUCEPTIBLE_CELL: {
            // If suceptible, count infected neighbors, decide if exposed
            size_t infectedNeighbors = get_infected_neighbors(b, x, y);
            int infectChance = (random() % 20) * infectedNeighbors;
            return (infectChance > 60) ? EXPOSED_CELL : SUCEPTIBLE_CELL;
	    }
        case EXPOSED_CELL:
            // Decide if moves to Infected
            // TODO Need to decide on a incubation period for the infection to start
            if (current_person->time_in_state >= 5) {
                return INFECTED_CELL;
            } else {}
            break;
	    case INFECTED_CELL:
            // Decide if moves to W or D or R
            // TODO Time will decide if => W, other will decide if D or R
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
    srandom(time(NULL) + world_rank);

    unsigned int num_rows = bc.height;
#if NUM_THREADS
	num_rows /= NUM_THREADS;
#endif
	unsigned int j = (*((size_t*)argp)) * num_rows;
	unsigned int end = j + num_rows - 1;
	for (; j < end; j++) {
		for (size_t i = 0; i < bc.width; i++) {
		    cell_state nextState = next_state(&bc, j, i);

		    if (bc.current[j][i].state != nextState) {
		        bc.next[j][i].time_in_state = 0;
		        bc.next[j][i].state = nextState;
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
    int rankUp = (world_rank == 0) ? (world_size - 1) : (world_rank - 1);
    int rankDown = (world_rank == (world_size - 1)) ? 0 : (world_rank + 1);

    // Copy Board
    for (size_t i = 0; i < b->height; i++) {
        memcpy(b->next[i], b->current[i], sizeof(Person*)*b->width);
    }

    send_above = b->current[0];
    send_below = b->current[b->height-1];

    // Send at the beginning of the tick
    // Send top of chunk to rank+1
    MPI_Isend(send_above, b->width, MPI_PERSON, rankUp, 0, MPI_COMM_WORLD, &sendUp);

    // Send bottom of chunk to rank-1
    MPI_Isend(send_below, b->width, MPI_PERSON, rankDown, 0, MPI_COMM_WORLD, &sendDown);

    // Recieve from rank-1 as bottom
    MPI_Irecv(recv_above, b->width, MPI_PERSON, rankUp, 0, MPI_COMM_WORLD, &recieveUp);

    // Recieve from rank+1 as top
    MPI_Irecv(recv_below, b->width, MPI_PERSON, rankDown, 0, MPI_COMM_WORLD, &recieveDown);


#if NUM_THREADS
	pthread_t* tid = calloc(NUM_THREADS, sizeof(pthread_t));
	size_t* tmpI = calloc(NUM_THREADS, sizeof(size_t));
#else
	size_t tmpI[] = {0};
#endif

    MPI_Wait(&sendUp, &mpi_stat);
    MPI_Wait(&sendDown, &mpi_stat);
    MPI_Wait(&recieveUp, &mpi_stat);
    MPI_Wait(&recieveDown, &mpi_stat);

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
	free(tmpI);
#endif

	Person** tmp = b->current;
	b->current = b->next;
	b->next = tmp;
}


int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    mpi_create_person_type();

    if (world_rank == 0) {
        g_start_cycles = (unsigned long long)GetTimeBase();
    }

	// Using built in random for now, may change out later
	srandom(time(NULL) + world_rank);

    size_t chunk = GRID_SIZE / world_size;
	InitBoard(&bc, GRID_SIZE, chunk);
    Send_Recv_Init(bc.width);

    unsigned int people, infected;

    people = count_people(&bc, SUCEPTIBLE_CELL);

    infect_people(&bc);

    infected = count_people(&bc, INFECTED_CELL);

    printf("Number of people suceptible: %d out of %d\n", people, GRID_SIZE*GRID_SIZE);
	printf("Number of people infected: %d out of %d\n", infected, people);

#ifdef DEBUG
    PrintBoard(&bc);
#endif

    unsigned long long infectedCount[MAX_TICKS];

    // Begin actual experiment
	for (size_t i = 0; i < MAX_TICKS; i++) {
	    tick(&bc);
        MPI_Barrier(MPI_COMM_WORLD);
        unsigned long long tmp = count_people(&bc, INFECTED_CELL) + count_people(&bc, WITHOUT_CELL);
        MPI_Reduce((void*)&tmp, (void*)&(infectedCount[i]), 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	}

    if (world_rank == 0) {
        for (size_t i = 0; i < MAX_TICKS; i++) {
            printf("Number of people infected: %lld out of %d\n", infectedCount[i], people);
        }
    }

    unsigned long long global_sum = 0;
    unsigned long long infect_no_spread = count_people(&bc, WITHOUT_CELL);
    MPI_Reduce((void*)&infect_no_spread, (void*)&global_sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    if (world_rank == 0) {
        infect_no_spread = global_sum;
        printf("Number of people infected without spreading: %lld out of %d\n", infect_no_spread, people);
    }

    global_sum = 0;
    unsigned long long recovered = count_people(&bc, RECOVERED_CELL);
    MPI_Reduce((void*)&recovered, (void*)&global_sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    if (world_rank == 0) {
        recovered = global_sum;
        printf("Number of people recovered: %lld out of %d\n", recovered, people);
    }

    global_sum = 0;
    unsigned long long dead = count_people(&bc, RECOVERED_CELL);
    MPI_Reduce((void*)&dead, (void*)&global_sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    if (world_rank == 0) {
        dead = global_sum;
        printf("Number of people dead: %lld out of %d\n", dead, people);
    }

#ifdef DEBUG
    PrintBoard(&bc);
#endif

    if (world_rank == 0) {
        g_end_cycles = (unsigned long long) GetTimeBase();
        g_time_in_secs = (double)(g_end_cycles - g_start_cycles) / PROC_FREQ;
        printf("Time = %f\n", g_time_in_secs);
    }

    DestroyBoard(&bc);
    //Send_Recv_Destroy();

	return 0;
}
