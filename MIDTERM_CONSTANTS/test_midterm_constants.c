#include "stdio.h"

#include "../MIDTERM_ERROR/midterm_error.h"
#include "midterm_constants.h"
#include "mpi.h"


int main(void){

	MPI_Init(0,0);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&num_proc);	
	initialize_grid_checker();

	printf("my_grid_col %d my_grid_row %d\n",my_grid_col,my_grid_row);	
	int rank = grid_index(my_grid_col, my_grid_row);

	printf("my_rank %d rank %d\n",my_rank,rank);
	MPI_Finalize();
	return 0;
}
