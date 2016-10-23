#include <stdio.h>

#include "mpi.h"

int main(int argc, char *argv[]){

	MPI_File fh;
	MPI_Datatype fileblk;
	MPI_Datatype ablk;
	MPI_Status status;
	int  rank;
	int num_proc;
	int NW = 5;
	int buf[NW*2];
	int i;
	for(i=0;i<(NW*2);i++){
		buf[i]=2147483647;
	}

	MPI_Init(0,0);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&num_proc);		

	MPI_File_open(MPI_COMM_WORLD,"test6.out",
	              MPI_MODE_CREATE | MPI_MODE_WRONLY,
                  MPI_INFO_NULL, &fh);

	/* Count, blocklen, stride */	
	MPI_Type_vector(2,NW,NW*num_proc,MPI_INT, &fileblk);
	MPI_Type_commit(&fileblk);

	MPI_Offset displacement = (MPI_Offset)rank*NW*sizeof(int);
	MPI_File_set_view(fh, displacement, MPI_INT, fileblk, "native", MPI_INFO_NULL);
	/* Count, type */
	MPI_Type_contiguous(NW, MPI_INT, &ablk);
	MPI_Type_commit(&ablk);
	MPI_File_write(fh, (void *)buf, 2, ablk, &status); 

	MPI_File_close(&fh);
	MPI_Type_free(&ablk);
	MPI_Type_free(&fileblk);
	MPI_Finalize();
	return 0;
}
