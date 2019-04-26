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

#define GRID_SIZE 4096
#define NUM_THREADS 16
// May look into making ticks represent hours, up to a number of days?
#define MAX_TICKS 1024
#define POPULATION_RATE 50 // Out of 100
#define INFECTION_RATE 1 // Out of 100
#define RECOVERY_RATE 3
#define INFECTIVE_TIME 8

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
    recv_above = (Person*)calloc(width, sizeof(Person));
    recv_below = (Person*)calloc(width, sizeof(Person));
}

void Send_Recv_Destroy() {
    // Each is one row of the grid
    free(recv_above);
    free(recv_below);
}

void InitPerson(Person* person) {
	person->age = (unsigned int)(random() % 89) + 1;
	person->time_in_state = 0;
    int chance = (int)(random() % 100);
    if (chance > POPULATION_RATE) {
        person->state = SUSCEPTIBLE_CELL;
    } else {
        person->state = FREE_CELL;
    }
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
            InitPerson(&(bc.current[i][j]));
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
		    //printf(" %d ", b->current[i][j].state);
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
			if (b->current[i][j].state != SUSCEPTIBLE_CELL) {
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
size_t get_infected_neighbors(Board* b, unsigned int y, unsigned int x) {
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
 * Decides if a given person will die from the disease based on their age.
 */
bool Compute_Death(Person* person){
    unsigned int age = person->age;
    srandom(time(NULL) + age);
    double death_chance = (random() % 1000) / 10.0;
    bool death = false;

    //Age range : 1 - 3, 8% mortality rate
    if (age > 1 && age <= 3) {
        if (death_chance < 8) {
            death = true;
        }
    }
    //Age range : 4 - 10, 5% mortality rate
    if (age > 4 && age <= 10) {
        if (death_chance < 5) {
            death = true;
        }
    }
    //Age range : 11 - 18, 2% mortality rate
    if (age > 11 && age <= 18) {
        if (death_chance < 2) {
            death = true;
        }
    }
    //Age range : 19 - 50, 0.5% mortality rate
    if (age > 19 && age <= 50) {
        if (death_chance < 0.5) {
            death = true;
        }
    }
    //Age range : 51 - 90, 2% mortality rate
    if (age > 51 && age <=60) {
        if (death_chance < 2) {
            death = true;
        }
    }
    //Age range : 61+ 8% mortality rate
    if (age > 51 && age <=60) {
        if (death_chance < 8) {
            death = true;
        }
    }

    return death;
}

/*
 * Computes the chance of recovery for a person who is infected
 * Uses the age, as well as the number of days they have been infected
 * If in state WITHOUT_CELL, add an additional 8 days (or however long the infective period is)
 */
bool Compute_Recovery(Person* person) {
    bool recovers = false;
    unsigned int age = person->age;
    unsigned int stateTime = person->time_in_state;
    srand(time(NULL) + (age * stateTime));
    if (person->state == WITHOUT_CELL) {
        stateTime += INFECTIVE_TIME;
    }

    if (stateTime > 15) {
        return true;
    } else {
        double recoveryChance = (random() % 1000) / 10.0;
        double recoveryRate = RECOVERY_RATE * stateTime;
        if (recoveryChance < recoveryRate) {
            recovers = true;
        }
    }

    return recovers;
}


int* getLiveNeighbors(Board* b, unsigned int x, unsigned int y) {
    int* neighbors = calloc(9, sizeof(int));
    int left = (x == 0) ? (b->width - 1) : (x - 1);
    int right = (x == b->width - 1) ? 0 : (x + 1);

    for (int i = 0; i < 9; i++) {
        neighbors[i] = 0;
    }

    /*
     * Return will be organized as follows
     *
     *      1  2  3
     *      4  5  6
     *      7  8  9
     *
     */
    // Only cells that cannot move are DEAD_CELL and FREE_CELL
    // All other are living and can move to a new state
    if (y == 0) {
        if (recv_below[x].state != DEAD_CELL && recv_below[x].state != FREE_CELL)// down 1 = 8
            neighbors[7] = 1;
        if (recv_below[left].state != DEAD_CELL && recv_below[left].state != FREE_CELL)
            neighbors[6] = 1; // Diagonal left down -> wrap to right edge
        if (recv_below[right].state != DEAD_CELL && recv_below[right].state != FREE_CELL)
           neighbors[8] = 1;// Diagonal right down
    }
    else {
        if (b->current[y-1][right].state != DEAD_CELL && b->current[y-1][right].state != FREE_CELL)
            neighbors[8] = 1;// down and right
        if (b->current[y-1][left].state != DEAD_CELL && b->current[y-1][left].state != FREE_CELL)
            neighbors[6] = 1;// down and left
        if (b->current[y-1][x].state != DEAD_CELL && b->current[y-1][x].state != FREE_CELL)
            neighbors[7] = 1;// down
    }


    if (y == b->height - 1) {
        if (recv_above[x].state != DEAD_CELL && recv_above[x].state != FREE_CELL)// up 1 = 2
            neighbors[1] = 1;
        if (recv_above[left].state != DEAD_CELL && recv_above[left].state != FREE_CELL)
            neighbors[0] = 1; // Diagonal left up -> wrap to right edge
        if (recv_above[right].state != DEAD_CELL && recv_above[right].state != FREE_CELL)
            neighbors[2] = 1;// Diagonal right up

    } else {
        if (b->current[y+1][right].state != DEAD_CELL && b->current[y+1][right].state != FREE_CELL)
            neighbors[2] = 1;// up and right
        if (b->current[y+1][left].state != DEAD_CELL && b->current[y+1][left].state != FREE_CELL)
            neighbors[0] = 1;// up and left
        if (b->current[y+1][x].state != DEAD_CELL && b->current[y+1][x].state != FREE_CELL)
            neighbors[1] = 1;// up

    }

    if (b->current[y][right].state != DEAD_CELL && b->current[y][right].state != FREE_CELL)
        neighbors[5] = 1;// just right  6
    if (b->current[y][left].state != DEAD_CELL && b->current[y][left].state != FREE_CELL)
        neighbors[3] = 1;// just left  4

    return neighbors;
}

/*
* This function will compute the next state for any given cell at the current time
* Parameters: x and y coordinates of target cell, and pointer to board itself
* Output: The value of the state that will succeed this state
*
*/
cell_state next_state(Board* b, unsigned int row, unsigned int col) {
    Person* current_person = &(b->current[row][col]);


	// Check current persons state to decide action
	switch (current_person->state) {
        case SUSCEPTIBLE_CELL: {
            // If suceptible, count infected neighbors, decide if exposed
            size_t infectedNeighbors = get_infected_neighbors(b, row, col);
            int infectChance = (random() % 20) * infectedNeighbors;
            return (infectChance > 15) ? EXPOSED_CELL : SUSCEPTIBLE_CELL;
        }
        case EXPOSED_CELL:
            // Decide if moves to Infected, wait time is 5 days
            if (current_person->time_in_state >= 5) {
                return INFECTED_CELL;
            } else {}
            break;
        case INFECTED_CELL:
            // Decide if moves to W or D or R
            if (current_person->time_in_state < INFECTIVE_TIME) {
                bool dies = Compute_Death(current_person);
                if (dies) {
                    return DEAD_CELL;
                }
                bool recovers = Compute_Recovery(current_person);
                if (recovers) {
                    return RECOVERED_CELL;
                }
                return INFECTED_CELL;
            }
            return WITHOUT_CELL;
        case WITHOUT_CELL: {

            // Decide if moves to R or D
            bool dies = Compute_Death(current_person);
            if (dies) {
                return DEAD_CELL;
            }
            bool recovers = Compute_Recovery(current_person);
            if (recovers) {
                return RECOVERED_CELL;
            }
            return WITHOUT_CELL;
        }
	    case FREE_CELL: {
            // Decide if person will move into free cell
            //int* people = getLiveNeighbors(b, x, y);
            //int move_chance = (random() % 9);
            //if (people[move_chance] == 1) {
              //  printf("would move a person here\n");
            //}
        }
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
	unsigned int row = (*((size_t*)argp)) * num_rows;
	unsigned int end = row + num_rows;
	for (; row < end; row++) {
		for (unsigned int col = 0; col < bc.width; col++) {
            cell_state nextState = next_state(&bc, row, col);
		    if (bc.current[row][col].state != nextState) {
		        bc.next[row][col].time_in_state = 0;
		        bc.next[row][col].state = nextState;
		    } else {
                bc.next[row][col].time_in_state++;
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
        memcpy(b->next[i], b->current[i], sizeof(Person) * b->width);
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

    unsigned long long people, infected, tmp;

    people = count_people(&bc, SUSCEPTIBLE_CELL);

    infect_people(&bc);

    infected = count_people(&bc, INFECTED_CELL);

    MPI_Reduce((void*)&people, (void*)&tmp, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    if (world_rank == 0) {
        people = tmp;
        //printf("Number of people susceptible: %lld out of %d\n", people, GRID_SIZE*GRID_SIZE);
    }

    MPI_Reduce((void*)&infected, (void*)&tmp, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    if (world_rank == 0) {
        infected = tmp;
        //printf("Number of people infected: %lld out of %lld\n", infected, people);
    }

    unsigned long long infectedCount[MAX_TICKS];
    //unsigned long long numInfected = (unsigned long long)infected;
    //unsigned long long recoveredCount[MAX_TICKS];

    // Begin actual experiment
    unsigned int day = 0;
    for (day = 0; day < MAX_TICKS; day++) {
#ifdef DEBUG
        PrintBoard(&bc);
#endif
	    tick(&bc);

        infectedCount[day] = count_people(&bc, INFECTED_CELL) + count_people(&bc, WITHOUT_CELL);

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce((void*)&(infectedCount[day]), (void*)&tmp, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        if (world_rank == 0) { infectedCount[day] = tmp; }
        if (tmp == 0) { break; }
        //tmp = count_people(&bc, RECOVERED_CELL);
        //MPI_Reduce((void*)&tmp, (void*)&(recoveredCount[i]), 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	}

    if (world_rank == 0) {
        printf("Days elapsed: %u\n", day);
        for (size_t i = 0; i < day; i++) {
            printf("%lld\n", infectedCount[i]);
            //printf("Number of people recovered: %lld out of %d\n", recoveredCount[i], people);
        }
    }

    unsigned long long global_sum = 0;
    unsigned long long recovered = count_people(&bc, RECOVERED_CELL);
    MPI_Reduce((void*)&recovered, (void*)&global_sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    if (world_rank == 0) {
        recovered = global_sum;
        //printf("Number of people recovered: %lld out of %lld\n", recovered, people);
    }

    global_sum = 0;
    unsigned long long dead = count_people(&bc, DEAD_CELL);
    MPI_Reduce((void*)&dead, (void*)&global_sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    if (world_rank == 0) {
        dead = global_sum;
       // printf("Number of people dead: %lld out of %lld\n", dead, people);
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
    Send_Recv_Destroy();

	return 0;
}
