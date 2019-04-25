//
// Created by cyrust on 4/25/2019.
//

#ifndef SEIR_BOARD_H
#define SEIR_BOARD_H

typedef struct {
    size_t width;
    size_t height;
    Person** current;
    Person** next;
} Board;

#endif //SEIR_BOARD_H
