#ifndef COMM_H
#define COMM_H 1

int comm_mpi_init(int* pargc, char*** pargv);
int comm_this_node(void);

#endif
