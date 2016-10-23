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

#include "mpi.h"

#include "../MIDTERM_CONSTANTS/midterm_constants.h"
#include "midterm_matrix.h"

int main(void){

	MPI_Init(0,0);

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

	if(num_proc!=2){
		fprintf(stderr,"ERROR only meant to test two processors\n");
		exit(1);
	}

	verbose = true;
	send_recv_sync = true;
	ghost_count = 1;
	matrix mat;
	initialize_grid_row_block();

	if(my_rank==0){
		mat = new_matrix_row_panel(4,5,0,1);
		print_matrix(mat);
		printf("Testing send_recv_ghost_rows\n");
		send_recv_ghost_rows(&mat);
		printf("After send_recv call\n");
		print_matrix(mat);	
	}else{
		mat = new_matrix_row_panel(4,5,1,0);
		set_all_elems(mat,1.0);
		send_recv_ghost_rows(&mat);
		print_matrix(mat);
	}

	delete_matrix(&mat);

	MPI_Finalize();	

	return 0;
}
