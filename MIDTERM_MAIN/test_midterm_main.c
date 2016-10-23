#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"

#include "../MIDTERM_MATRIX/midterm_matrix.h"
#include "../MIDTERM_IO/midterm_io.h"
#include "../MIDTERM_CONSTANTS/midterm_constants.h"
#include "../MIDTERM_ERROR/midterm_error.h"

int main(int argc, char **argv){

	float local_bug_sum;
	int counter;
	int i;
	int j;
	float val;

	#ifdef _MPI_H_
	MPI_Init(0,0);
	
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	
	MPI_Comm_size(MPI_COMM_WORLD,&num_proc);

	options(argc, argv);

	initialize_grid_row_block();

	//create_filename_out(my_rank);
	//switch_stdout(filename_out);
	matrix mat_a = NULL;
	matrix mat_b = NULL;
	matrix * temp = (matrix *) malloc(sizeof(matrix *));

	//read_pgm(&mat_a, &mat_b);
	if (my_rank==0 && num_proc==1){
		mat_a = new_matrix(4,4);
		mat_b = new_matrix(4,4);
		set_all_elems(mat_a,1.0);
/*
		set_elem_core_matrix(mat_a,0,0,1.0);
		set_elem_core_matrix(mat_a,3,0,1.0);
		set_elem_core_matrix(mat_a,2,0,1.0);
		set_elem_core_matrix(mat_a,2,1,1.0);
		set_elem_core_matrix(mat_a,2,3,1.0);
*/
	} else {
		if(my_rank==0){
			mat_a = new_matrix_row_panel( 2,4,0,1);
			mat_b = new_matrix_row_panel( 2,4,0,1);
			set_all_elems(mat_a,1.0);
//			set_elem_core_matrix(mat_a,0,0,1.0);
		}else{
			mat_a = new_matrix_row_panel( 2,4,1,0);
			mat_b = new_matrix_row_panel( 2,4,1,0);
			set_all_elems(mat_a,0.0);
/*			set_elem_core_matrix(mat_a,1,0,1.0);
			set_elem_core_matrix(mat_a,0,0,1.0);
			set_elem_core_matrix(mat_a,0,1,1.0);
			set_elem_core_matrix(mat_a,0,3,1.0);
*/
		}
	}


	*temp = mat_a;

	//print_matrix(*temp);

	free(init_file);

	/* Initialize the number of bugs in the system */
	global_bug_sum = 0;	
	
	/* Begin Conway's Game of Life */
	for( counter=0; counter<init_iter; counter++){
	

		if(counter==0){
			send_recv_ghost_rows(temp);
			print_matrix(*temp);
		}	
		local_bug_sum = sum_all_core_matrix_elems(*temp);
		printf("Local sum %f\n",local_bug_sum);
		MPI_Reduce(&local_bug_sum,&global_bug_sum,1,MPI_FLOAT,MPI_SUM, 0,MPI_COMM_WORLD);
		if(my_rank==0){
			printf("Global sum %f\n",global_bug_sum);
		}
		/* Alternative which matrix is considered new */
		if(counter%2==1){
			//print_matrix(mat_a);
			//print_matrix(mat_b);
			exit(1);
			//conways_game_block(temp, &mat_a);
			*temp = mat_a;
		}else{
			//conways_game_block(temp, &mat_b);
			*temp = mat_b;	
		}

	}

	delete_matrix(&mat_a);
	delete_matrix(&mat_b);
	free(temp);
	//revert_stdout();
	MPI_Finalize();	
	#endif
	return 0;
}
