#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "../MIDTERM_CONSTANTS/midterm_constants.h"
#include "../MIDTERM_ERROR/midterm_error.h"
#include "../MIDTERM_MATRIX/midterm_matrix.h"

#include "mpi.h"
int main(int argc, char **argv){

	MPI_Init(0,0);
	
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	
	MPI_Comm_size(MPI_COMM_WORLD,&num_proc);

	options(argc, argv);

	if(block_on){
		initialize_grid_row_block();
	}else{
		initialize_grid_checker();
	}
	matrix mat_a = NULL;
	matrix mat_b = NULL;	

	printf("init_file in main %s\n",init_file);

	read_pgm(&mat_a,&mat_b);	

	write_pgm(mat_a, 1, my_rank);

	free(init_file);

	delete_matrix(&mat_a);
	delete_matrix(&mat_b);


	MPI_Finalize();

	printf("Testing String functions\n");

//	char * str1 = "Hello";
//	int rv;
//	rv = string_check_for_dashs(str1);
//	assert(rv==0);
//	char * str2 = "Hi,Yo";
//	rv = string_check_for_dashs(str2);
	//assert(rv==1);
	
	return 0;
}
