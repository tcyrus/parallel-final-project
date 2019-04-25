//
// Created by cyrust on 4/25/2019.
//

#include <mpi.h>
#include <stddef.h>

#include "cell_state.h"

#ifndef SEIR_PERSON_H
#define SEIR_PERSON_H

typedef struct {
    unsigned int time_in_state;
    cell_state state;
    unsigned int age;
} Person;

MPI_Datatype MPI_PERSON;

void mpi_create_person_type() {
    const int nitems = 3;
    int blocklengths[] = {1, 1, 1};
    MPI_Datatype types[] = {MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED};
    MPI_Aint offsets[] = {
            (MPI_Aint)offsetof(Person, time_in_state),
            (MPI_Aint)offsetof(Person, state),
            (MPI_Aint)offsetof(Person, age)
    };
    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &MPI_PERSON);
    MPI_Type_commit(&MPI_PERSON);
}

#endif //SEIR_PERSON_H
