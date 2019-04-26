//
// Created by cyrust on 4/25/2019.
//

#ifndef SEIR_CELL_STATE_H
#define SEIR_CELL_STATE_H

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
/*
typedef enum {
    BOARDER_CELL,
    FREE_CELL,
    SUSCEPTIBLE_CELL,
    EXPOSED_CELL,
    INFECTED_CELL,
    RECOVERED_CELL,
    DEAD_CELL,
    WITHOUT_CELL
} cell_state;
*/

#define cell_state unsigned int
#define BOARDER_CELL 0
#define FREE_CELL 1
#define SUSCEPTIBLE_CELL 2
#define EXPOSED_CELL 3
#define INFECTED_CELL 4
#define RECOVERED_CELL 5
#define DEAD_CELL 6
#define WITHOUT_CELL 7

#ifdef DEBUG
// In case we actually want to print out the letters instead of numbers
const char cell_state_names[] = { 'B', 'F', 'S', 'E', 'I', 'R', 'D', 'W' };
#endif


#endif //SEIR_CELL_STATE_H
