#include <mpi.h>
#include <fstream>
#include <iostream>
using namespace std;

int main(int argc, char const *argv[])
{
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	int M, N;
	int rank, numprocs, ierr;

	
}