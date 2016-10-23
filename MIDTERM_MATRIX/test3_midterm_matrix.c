/*
 *          Midterm 3
 *          Joshua Brown
 *          Scientific Computing
 *          test2_midterm_matrix.c function
 *          Conway's Game of Life
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "../MIDTERM_ERROR/midterm_error.h"
#ifdef _MPI_H_
#include "mpi.h"
#endif
#include "../MIDTERM_CONSTANTS/midterm_constants.h"
#include "midterm_matrix.h"

int main(void){

	MPI_Init(0,0);

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

	if(num_proc!=3){
		fprintf(stderr,"ERROR only meant to test three processors\n");
		exit(1);
	}


	verbose = true;
	send_recv_sync = true;
	ghost_count = 1;
	initialize_grid_row_block();
	matrix mat;

	float local_sum;
	float global_sum;

	create_filename_out(my_rank);
	switch_stdout(filename_out);

	if(my_rank==0){
		mat = new_matrix_row_panel(4,5,0,1);
		printf("Testing send_recv_ghost_rows\n");
		send_recv_ghost_rows(&mat);
		printf("After send_recv call\n");
		print_matrix(mat);	
	}else if(my_rank==(num_proc-1)){
		mat = new_matrix_row_panel(4,5,1,0);
		set_all_elems(mat,2.0);
		send_recv_ghost_rows(&mat);
		print_matrix(mat);
	}else{
		mat = new_matrix_row_panel(4,5,1,1);
		set_all_elems(mat,1.0);
		send_recv_ghost_rows(&mat);
		print_matrix(mat);
	}

	local_sum = sum_all_core_matrix_elems(mat);
	MPI_Reduce(&local_sum, &global_sum, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

	if(my_rank==0){
		printf("global sum %f \n",global_sum);
	}
	delete_matrix(&mat);

	revert_stdout();
	MPI_Finalize();	

	return 0;
}
