/*
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
#include "../MIDTERM_ERROR/midterm_error.h"
#include "midterm_matrix.h"

int main(void){

	MPI_Init(0,0);

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

	if(num_proc!=4){
		fprintf(stderr,"ERROR only meant to test four processors\n");
		exit(1);
	}

	verbose = true;
	send_recv_sync=true;
	ghost_count = 1;
	matrix mat;

	float local_sum;
	float global_sum;

	create_filename_out(my_rank);
	switch_stdout(filename_out);

	initialize_grid_checker();

	printf("Creating Matrx\n");
	printf("my_grid_row %d my_grid_col %d\n",my_grid_row,my_grid_col);
	printf("grid_dim_row %d grid_dim_col %d\n",grid_dim_row,grid_dim_col);
	printf("My rank from grid_index %d\n",grid_index(my_grid_row,my_grid_col));
	if(my_rank==0){
		mat = new_matrix_checkerboard_elem(4,4,0,1,1,0);
	}else if(my_rank==1){
		mat = new_matrix_checkerboard_elem(4,4,0,1,0,1);
		set_all_elems(mat,1.0);
	}else if(my_rank==2){
		mat = new_matrix_checkerboard_elem(4,4,1,0,1,0);
		set_all_elems(mat,2.0);
	}else{
		mat = new_matrix_checkerboard_elem(4,4,1,0,0,1);
		set_all_elems(mat,3.0);
	}

	print_matrix(mat);	
	printf("Testing send_recv_ghost_cols\n");
	send_recv_ghost_cols(&mat);
	print_matrix(mat);	
	printf("Testing send_recv_ghost_rows\n");
	send_recv_ghost_rows(&mat);
	print_matrix(mat);
	local_sum = sum_all_core_matrix_elems(mat);
	printf("local sum %f\n",local_sum);
	MPI_Reduce(&local_sum, &global_sum, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

	if(my_rank==0){
		printf("global sum %f \n",global_sum);
	}
	delete_matrix(&mat);

	revert_stdout();
	MPI_Finalize();	

	return 0;
}
