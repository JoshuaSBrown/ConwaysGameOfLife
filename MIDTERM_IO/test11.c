#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "mpi.h"

int main(void){

	int rank;
	int num_proc;
	int i;
	
	int NW = 3;
	float array[NW*2];
	char char_array[NW*2];

	MPI_Init(0,0);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&num_proc);

	if(rank==0){
		array[0] = 0.0;
		array[1] = 0.0;
		array[2] = 0.0;
		array[3] = 0.0;
		array[4] = 0.0;
		array[5] = 0.0;
	}else{
		array[0] = 1.0;
		array[1] = 1.0;
		array[2] = 1.0;
		array[3] = 1.0;
		array[4] = 1.0;
		array[5] = 1.0; 
	}

	for(i=0;i<(NW*2);i++){
		char_array[i]=((int)array[i])*255;
	}
	MPI_File fh;
    MPI_Datatype fileblk;
	MPI_Datatype ablk;
	MPI_Datatype float_array;
	MPI_Status status;

	MPI_File_open(MPI_COMM_WORLD,"test11.out",
                  MPI_MODE_CREATE | MPI_MODE_WRONLY,
                  MPI_INFO_NULL, &fh);
	
	MPI_Type_vector(2,NW,NW*num_proc,MPI_CHAR,&fileblk);
	MPI_Type_commit(&fileblk);

	MPI_Type_contiguous(NW, MPI_CHAR, &ablk);
	MPI_Type_commit(&ablk);

	MPI_Offset displacement = (MPI_Offset)(rank*NW*sizeof(char));
	MPI_File_set_view(fh,displacement,MPI_INT,fileblk,"native",MPI_INFO_NULL);

	MPI_File_write(fh,char_array,2,ablk,&status);

	MPI_File_close(&fh);
	MPI_Type_free(&fileblk);
	MPI_Type_free(&ablk);
	MPI_Finalize();
	return 0;
}
