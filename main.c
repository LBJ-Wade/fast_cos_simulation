#include <stdio.h>
#include "comm.h"
#include "msg.h"

int main(int argc, char* argv[])
{
  // Setup MPI Init  as comm_mpi_init?
  //
  comm_mpi_init(&argc, &argv);

  msg_printf(msg_info, "Hello World\n");

  comm_mpi_finalise();
}
