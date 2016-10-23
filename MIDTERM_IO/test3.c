#include <stdio.h>
#include "mpi.h"

MPI_Status status;
MPI_File fh;
MPI_Offset offset;
int rank;

int main(int argc, char * argv){

	MPI_Init(0,0);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	MPI_File_open(MPI_COMM_WORLD, "sample.txt",
				  MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

	nints = FILESIZE / (nprocs*INTSIZE);

	MPI_File_read_at(fh, offset, buf, nints, MPI_INT, &status);

	MPI_Get_count(&status, MPI_INT, &count);

	MPI_File_close(&fh);

	MPI_Finalize();
	return 0;
}
