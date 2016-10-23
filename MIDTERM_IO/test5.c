#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

int main(int argc, char *argv[]){


	MPI_File fh;
	MPI_Datatype fileblk;
	MPI_Datatype ablk;
	MPI_Status status;
	int rank;
	int num_proc;
	int NW = 10;
	char buf2[NW*2];
	int i;

	MPI_Init(0,0);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&num_proc);
	
	if(rank==0){	 
		sprintf(buf2, "First Row\nThird Row\n"); 

	}else{
		sprintf(buf2, "SecondRow\nFourthRow\n"); 
	}

	MPI_File_open(MPI_COMM_WORLD,"test5.out",
                  MPI_MODE_CREATE | MPI_MODE_WRONLY,
                  MPI_INFO_NULL, &fh);

	/* Define file type */ 
    MPI_Type_vector(2,NW,NW*num_proc,MPI_CHAR, &fileblk);
	MPI_Type_commit(&fileblk);

	/* Define file view */
	MPI_Offset displacement = (MPI_Offset)rank*NW*sizeof(char);
	MPI_File_set_view(fh,displacement,MPI_INT,fileblk,"native",MPI_INFO_NULL);

	/* Define e type */
	MPI_Type_contiguous(NW, MPI_CHAR, &ablk);
	MPI_Type_commit(&ablk);
	
	/* Write rows to the file */
	MPI_File_write(fh,(void *)buf2, 2, ablk, &status);

	MPI_File_close(&fh);
	MPI_Type_free(&ablk);
	MPI_Type_free(&fileblk);
	MPI_Finalize();
	return 0;
}
