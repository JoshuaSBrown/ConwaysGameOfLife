#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpi.h"

int main(void){

	MPI_Status status;
	MPI_File fh;
	MPI_Datatype filetype;
	int my_rank;
	int num_proc;
	int ndims;
	int array_of_gsizes[2];
	int array_of_distribs[2];
	int array_of_dargs[2];
	int array_of_psizes[2];

	ndims = 2;
	
	int NW = 10;
	char buf2[NW*2];

	MPI_Init(0,0);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&num_proc);

	if(my_rank==0){
		sprintf(buf2, "First Row\nThird Row\n");

	}else{
		sprintf(buf2, "SecondRow\nFourthRow\n");
	}

	/* Number of rows columns etc 
	 * How large the dimensions are */
	array_of_gsizes[0] = 40;
	array_of_gsizes[1] = 1;
	
	/* Don't know */
	array_of_distribs[0] = MPI_DISTRIBUTE_BLOCK;
	array_of_distribs[1] = MPI_DISTRIBUTE_BLOCK;
	
	/* Don't know */
	array_of_dargs[0] = 20;
	array_of_dargs[1] = MPI_DISTRIBUTE_DFLT_DARG;

	/* Number of processors in each dimension */
	array_of_psizes[0] = 2;
	array_of_psizes[1] = 1;
	
	MPI_Type_create_darray(num_proc,
                           my_rank,
                           ndims,
                           array_of_gsizes,
                           array_of_distribs,
                           array_of_dargs,
                           array_of_psizes,
                           MPI_ORDER_C,
                           MPI_CHAR, 
                           &filetype);

	MPI_Type_commit(&filetype);

	MPI_File_open(MPI_COMM_WORLD,"sample12.txt",
                  MPI_MODE_CREATE | MPI_MODE_WRONLY,
                  MPI_INFO_NULL, &fh);
	
	MPI_File_set_view(fh, 0, MPI_CHAR, filetype, "native", 
                      MPI_INFO_NULL);

	MPI_File_write_all(fh,buf2,strlen(buf2),MPI_CHAR,&status);
	MPI_File_close(&fh);

	MPI_Finalize();

	return 0;
}
