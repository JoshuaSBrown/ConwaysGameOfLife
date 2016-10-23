#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "mpi.h"

int main(void){

	int rank;
	int num_proc;


    int sizeoffloat = (int)sizeof(float);
    int sizeofchar  = (int)sizeof(char);
    int stride_conversion = sizeoffloat/sizeofchar;

	printf("stride %d\n",stride_conversion);
	
	FILE * fp;

	float val = 1.0;
	fp = fopen("test9.out","w");
	
	fputc(val,fp);
	close(fp);	

	int NW = 3;
	float array[NW*2];

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
		array[0] = FLT_MAX;
		array[1] = FLT_MAX;
		array[2] = FLT_MAX;
		array[3] = FLT_MAX;
		array[4] = FLT_MAX;
		array[5] = FLT_MAX;
	}

	MPI_File fh;
    MPI_Datatype fileblk;
	MPI_Datatype ablk;
	MPI_Datatype float_array;
	MPI_Status status;

	MPI_File_open(MPI_COMM_WORLD,"test9b.out",
                  MPI_MODE_CREATE | MPI_MODE_WRONLY,
                  MPI_INFO_NULL, &fh);
	/* Count 3 of them block lenth 1 stride conversion should be 4
     * Defining the etype*/
	//MPI_Type_vector(3,1,stride_conversion, MPI_CHAR,&ablk);
	//MPI_Type_commit(&ablk);

	MPI_Type_contiguous(6,MPI_FLOAT,&float_array);
	MPI_Type_commit(&float_array);
	//MPI_Type_vector(2,1,num_proc,ablk, &fileblk);
	//MPI_Type_commit(&fileblk);
	/* Define file view */
	//MPI_Offset displacement = (MPI_Offset)(rank*NW*(sizeof(char)));
	MPI_File_set_view(fh,0,MPI_FLOAT,float_array,"native",MPI_INFO_NULL);

	MPI_File_write(fh,(void *)array,6,MPI_FLOAT,&status);

	MPI_File_close(&fh);
	//MPI_Type_free(&fileblk);
	//MPI_Type_free(&ablk);
	MPI_Finalize();
	return 0;
}
