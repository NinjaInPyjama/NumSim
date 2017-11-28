#include <stdio.h>
//#include <mpich/mpi.h>
#include <mpi.h>

int main (int argc, char **argv) {
  int _rank = 0;
  int _size = 0;
  
  // Initialize MPI and pass all parameters down to all threads
  MPI_Init(&argc, &argv);
  
  // Get the rank of the current thread
  MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
  // Get total number of threads
	MPI_Comm_size(MPI_COMM_WORLD, &_size);
  
  // Print a stupid message. Caution: Lot's of threads might cause a resource conflict
  printf("There %s %i World%s. Hello World %i!\n", (_size==1?"is":"are"), _size, (_size==1?"":"s"), _rank);
  
  // Wait for everyone
  MPI_Barrier(MPI_COMM_WORLD);
  
  // I'm here to store information nobody wants to know -.-
  MPI_Status stat;
  
  // Some space to talk to each other
  double buffer[42];
  for (int i=0;i<42;++i) buffer[i] = _rank*42 + i; // Such unique, much wow!
  
  printf("Hello World %i, look at this: %.2lf %.2lf\n", _rank, buffer[0], buffer[41]);
  
  // Wait for everyone
  MPI_Barrier(MPI_COMM_WORLD);
  
  // To, From
  int dest = _rank ^ 1; // So 0 sends to 1, 1 to 0, 2 -> 3, 3 -> 2, ...
  // This is our communication tag. It can be anything as long both parts have the same tag
  int tag = 4;
  
  // No send-receive the s**t out of it
  MPI_Sendrecv_replace(buffer, 42, MPI_DOUBLE, dest, tag, dest, tag, MPI_COMM_WORLD, &stat);
  // send and receive from/to buffer
  // send/receive 42 elements
  // the elements are MPI_DOUBLE aka double, thus 8 byte
  // send to dest
  // with tag tag
  // receive from dest
  // with tag tag
  // don't care abaout the rest. It has been there since ever, and it will stay.
  
  printf("Hello World %i, look at this: %.2lf %.2lf\n", _rank, buffer[0], buffer[41]);
  
  // Wait for everyone
  MPI_Barrier(MPI_COMM_WORLD);
  
  // And now to all the even threads: Put your hands up in the air!
  if (_rank%2 == 0) {
    MPI_Send(buffer, 42, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
  } else { // And to all the odd threads: Look at those hands!
    MPI_Recv(buffer, 42, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, &stat);
  }
  
  printf("Hello World %i, look at this: %.2lf %.2lf\n", _rank, buffer[0], buffer[41]);
  
  // Wait for everyone
  MPI_Barrier(MPI_COMM_WORLD);
  
  // And now to something completly different
  double myvalue = _rank;
  double ourvalue;
  
  MPI_Allreduce(&myvalue, &ourvalue, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  // Send myvalue
  // Write to ourvalue
  // the array is 1 element long
  // We are still using doubles
  // Perform a SUM with all the data provided by all threads. Some prefere MPI_MIN or MPI_MAX
  // Leave the MPI_COMM_WORLD alone. It dislikes people
  
  printf("Hello World %i, look at this: %.2lf\n", _rank, ourvalue);
  
  // Wait for everyone
  MPI_Barrier(MPI_COMM_WORLD);
  
  // Congratulation! You're done!
  MPI_Finalize();
  return 0;
}
