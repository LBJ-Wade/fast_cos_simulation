#ifndef COMM_H
#define COMM_H 1

void comm_mpi_init(int* pargc, char*** pargv);
void comm_mpi_finalise(void);
void comm_abort(void);
int comm_this_node(void);
int comm_n_nodes(void);
void comm_bcast_int(int* p_int, int count);
void comm_bcast_double(double* p_double, int count);
#endif
