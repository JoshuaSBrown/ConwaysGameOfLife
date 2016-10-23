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
#include "../MIDTERM_IO/midterm_io.h"
#include "../MIDTERM_ERROR/midterm_error.h"
#include "midterm_matrix.h"

int main(void){

	MPI_Init(0,0);

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

	initialize_grid_checker();

	init_test = 1;
	verbose = true;
	send_recv_sync=true;
	ghost_count = 1;
	matrix mat_a;
	matrix mat_b;
	matrix extra;

	float local_sum;
	float global_sum;
	int width;
	int height;

	width  = grid_dim_col*12;
	height = grid_dim_row*12;
	sub_matrix_dim_columns = width  / grid_dim_col;
	sub_matrix_dim_rows    = height / grid_dim_row;
	printf("sub_col  %d\n",sub_matrix_dim_columns);
	printf("width    %d\n",width);
	printf("grid_col %d\n",grid_dim_col);
	printf("sub_row  %d\n",sub_matrix_dim_rows);
	printf("height   %d\n",height);
	printf("grid_row %d\n",grid_dim_row);

	mat_a = define_matrix_blockchecker();
	mat_b = define_matrix_blockchecker();
	extra = define_matrix_blockchecker();

	set_all_elems(mat_a,my_rank);
	set_all_elems(mat_b,my_rank);

	create_filename_out(my_rank);
	switch_stdout(filename_out);

	printf("Matrices\n");
	printf("*****************************************\n");
	fflush(NULL);
	if(num_proc>1){
		send_recv_ghost_cols(&mat_a);
		print_matrix_compressed(mat_a);
		send_recv_ghost_rows(&mat_a);
		print_matrix_compressed(mat_a);
	}
	printf("Matrices\n");
	printf("*****************************************\n");
	fflush(NULL);
	conways_game(&mat_a,&mat_b,&extra);

	print_matrix_compressed(extra);

	if(num_proc>1){
		send_recv_ghost_cols(&extra);
		print_matrix_compressed(extra);
		send_recv_ghost_rows(&extra);
		print_matrix_compressed(extra);
	}
	conways_game(&extra,&mat_b,&mat_a);
	print_matrix_compressed(mat_a);
	revert_stdout();
	MPI_Finalize();	

	return 0;
}
