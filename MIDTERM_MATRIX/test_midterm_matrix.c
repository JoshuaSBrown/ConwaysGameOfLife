/*
 *          Joshua Brown
 *          Scientific Computing
 *          test_midterm_matrix.c function
 *          Conway's Game of Life
 *
 * This function tests the following functions:
 *		o new_matrix
 *      o print_matrix
 *      o delete_matrix
 *      o get_rows
 *      o get_cols			
 *		o get_total_elems
 *		o new_matrix_row_panel
 *      o set_elem_core_matrix
 *      o get_elem_core_matrix
 *      o set_elem_top_row_ghost_matrix
 *		o get_elem_top_row_ghost_matrix      
 *      o set_elem_bottom_row_ghost_matrix
 *		o get_elem_bottom_row_ghost_matrix      
 *      o get_rows_core_matrix
 *      o get_cols_core_matrix
 */
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "../MIDTERM_CONSTANTS/midterm_constants.h"
#include "../MIDTERM_ERROR/midterm_error.h"
#ifdef _MPI_H_
#include "mpi.h"
#endif
#include "midterm_matrix.h"


int main(void){

	int rv;
	float rvd;
	bool rvbool;
	matrix mat;
	matrix new_mat;
	matrix checker;
	int i;
	int j;

	printf("Testing: new_matrix\n");
	mat = new_matrix(0,1);
	assert(mat==NULL);
	mat = new_matrix(1,0); assert(mat==NULL);
	mat = new_matrix(1,1);
	assert(mat!=NULL);

	printf("Testing: print_matrix\n");
	rv = print_matrix(NULL);
	assert(rv==-1);
	rv = print_matrix(mat);
	assert(rv==0);

	printf("Testing: delete_matrix\n");
	rv = delete_matrix(NULL);
	assert(rv==-1);
	rv = delete_matrix(&mat);
	assert(rv==0);
	
	mat = new_matrix(2,3);
	assert(mat!=NULL);
	printf("Testing: get_rows\n");
	rv = get_rows(NULL);
	assert(rv==-1);
	rv = get_rows(mat);
	assert(rv==2);
	
	printf("Testing: get_cols\n");
	rv = get_cols(NULL);
	assert(rv==-1);
	rv = get_cols(mat);
	assert(rv==3);
	
	printf("Testing: get_total_elems\n");
	rv = get_total_elems(NULL);
	assert(rv==-1);
	rv = get_total_elems(mat);
	assert(rv==(2*3));
	delete_matrix(&mat);

	printf("Testing: new_matrix_row_panel\n");
	mat = new_matrix_row_panel(0, 1, 0, 0);
	assert(mat==NULL);
	mat = new_matrix_row_panel(1, 0, 0, 0);
	assert(mat==NULL);
	mat = new_matrix_row_panel(1, 1, -1, 0);
	assert(mat==NULL);
	mat = new_matrix_row_panel(1, 1, 0, -1);
	assert(mat==NULL);
	mat = new_matrix_row_panel(1, 1, 0, 0);
	assert(mat!=NULL);
	printf("Only Core matrix 1,1\n");
	print_matrix(mat);
	rv = delete_matrix(&mat);
	assert(rv == 0);
	assert(mat==NULL);
	printf("Ghost row on the top of a 2,2 matrix\n");
	mat = new_matrix_row_panel(2, 2, 1, 0);
	assert(mat!=NULL);
	rv = print_matrix(mat);
	assert(rv==0);
	rv = delete_matrix(&mat);
	assert(rv==0);
	assert(mat==NULL);
	printf("Ghost row on the bottom of a 3,3 matrix\n");
	mat = new_matrix_row_panel(3, 3, 0, 1);
	assert(mat!=NULL);
	rv = print_matrix(mat);
	assert(rv == 0);
	rv = delete_matrix(&mat);
	assert(rv==0);
	assert(mat==NULL);
	printf("Ghost row above and below a 2,2 matrix\n");
	mat = new_matrix_row_panel(2, 2, 1, 1);
	assert(mat!=NULL);
	rv = print_matrix(mat);
	assert(rv==0);
	rv = delete_matrix(&mat);
	assert(rv==0);
	assert(mat==NULL);

	printf("Testing: set_elem_core_matrix & get_elem_core_matrix\n");
	mat = new_matrix_row_panel(3, 2, 1, 1);
	assert(mat!=NULL);
	print_matrix(mat);
	rv = set_elem_core_matrix(NULL,0,0,2.0);
	assert(rv==-1);
	rv = set_elem_core_matrix(mat, -1,0,2.0);
	assert(rv==-1);
	rv = set_elem_core_matrix(mat, 0,-1,2.0);
	assert(rv==-1);
	rv = set_elem_core_matrix(mat, 0, 0,2.0);
	assert(rv==0);
	
	rvd = get_elem_core_matrix(mat, 0, 0);
	assert(rvd==2.0);

	rv = set_elem_core_matrix(mat, 3, 1, 3.0);	
	assert(rv==-1);
	rv = set_elem_core_matrix(mat, 2, 2, 5.0);
	assert(rv==-1);
	rv = set_elem_core_matrix(mat, 2, 1, 8.0);
	assert(rv==0);
	rvd = get_elem_core_matrix(mat, 2, 1);
	assert(rvd==8.0);

	rvd = get_elem_core_matrix(NULL, 2, 1);
	assert(rvd==-1.0);
	rvd = get_elem_core_matrix(mat, 3,1);
	assert(rvd==-1.0);
	rvd = get_elem_core_matrix(mat, 2,2);
	assert(rvd==-1.0);
	rvd = get_elem_core_matrix(mat, -1, 0);
	assert(rvd==-1.0);
	rvd = get_elem_core_matrix(mat, 0, -1);
	assert(rvd==-1.0);
	rvd = get_elem_core_matrix(mat, 0, 0);
	assert(rvd==2.0);
	print_matrix(mat);	

	printf("Testing: set_elem_top_row_ghost_matrix\n");
	rv = set_elem_top_row_ghost_matrix(NULL,0,0,21.3);
	assert(rv==-1);
	rv = set_elem_top_row_ghost_matrix(mat,1,0,23.8);
	assert(rv==-1);
	rv = set_elem_top_row_ghost_matrix(mat,0,2,98);
	assert(rv==-1);
	rv = set_elem_top_row_ghost_matrix(mat,0,1,7);
	assert(rv==0);
	rv = set_elem_top_row_ghost_matrix(mat,-1,0,54);
	assert(rv==-1);
	rv = set_elem_top_row_ghost_matrix(mat,0,-1,42);
	assert(rv==-1);
	rv = set_elem_top_row_ghost_matrix(mat,0,0,84);
	assert(rv==0);

	printf("Testing: get_elem_top_row_ghost_matrix\n");
	rvd = get_elem_top_row_ghost_matrix(NULL,0,0);
	assert(rvd==-1.0);
	rvd = get_elem_top_row_ghost_matrix(mat,-1,0);
	assert(rvd==-1.0);
	rvd = get_elem_top_row_ghost_matrix(mat,0,-1);
	assert(rvd==-1.0);
	rvd = get_elem_top_row_ghost_matrix(mat,1,0);
	assert(rvd==-1);
	rvd = get_elem_top_row_ghost_matrix(mat,0,2);
	assert(rvd==-1.0);
	rvd = get_elem_top_row_ghost_matrix(mat,0,1);
	assert(rvd==7);
	rvd = get_elem_top_row_ghost_matrix(mat,0,0);
	assert(rvd==84);

	printf("Testing: set_elem_bottom_row_ghost_matrix\n");
	rv = set_elem_bottom_row_ghost_matrix(NULL,0,0,21.3);
	printf("1\n");
	assert(rv==-1);
	rv = set_elem_bottom_row_ghost_matrix(mat,1,0,23.8);
	printf("2\n");
	assert(rv==-1);
	rv = set_elem_bottom_row_ghost_matrix(mat,0,2,98);
	printf("3\n");
	assert(rv==-1);
	rv = set_elem_bottom_row_ghost_matrix(mat,0,0,7);
	printf("4\n");
	assert(rv==0);
	rv = set_elem_bottom_row_ghost_matrix(mat,-1,0,54);
	printf("5\n");
	assert(rv==-1);
	rv = set_elem_bottom_row_ghost_matrix(mat,0,-1,42);
	printf("6\n");
	assert(rv==-1);
	rv = set_elem_bottom_row_ghost_matrix(mat,0,1,84);
	printf("7\n");
	assert(rv==0);

	printf("Testing: get_elem_bottom_row_ghost_matrix\n");
	rvd = get_elem_bottom_row_ghost_matrix(NULL,0,0);
	assert(rvd==-1.0);
	rvd = get_elem_bottom_row_ghost_matrix(mat,-1,0);
	assert(rvd==-1.0);
	rvd = get_elem_bottom_row_ghost_matrix(mat,0,-1);
	assert(rvd==-1.0);
	rvd = get_elem_bottom_row_ghost_matrix(mat,1,0);
	assert(rvd==-1);
	rvd = get_elem_bottom_row_ghost_matrix(mat,0,2);
	assert(rvd==-1.0);
	rvd = get_elem_bottom_row_ghost_matrix(mat,0,0);
	assert(rvd==7);
	rvd = get_elem_bottom_row_ghost_matrix(mat,0,1);
	assert(rvd==84);
	print_matrix(mat);

	printf("Testing: get_rows_core_matrix\n");
	rv = get_rows_core_matrix(NULL);
	assert(rv==-1);
	rv = get_rows_core_matrix(mat);
	assert(rv==3);
	
	printf("Testing: get_cols_core_matrix\n");
	rv = get_cols_core_matrix(NULL);
	assert(rv==-1);
	rv = get_cols_core_matrix(mat);
	assert(rv==2);

	printf("Testing: sum_all_core_matrix_elems\n");
	rvd = sum_all_core_matrix_elems(NULL);
	assert(rvd==-1.0);
	rvd = sum_all_core_matrix_elems(mat);
	printf("sum %f\n",rvd);
	assert((int)rvd==10);
	
	printf("Testing: within_core_matrix\n");
	print_matrix(mat);
	rvbool = within_core_matrix(NULL,1,1);
	assert(rvbool==false);
	rvbool = within_core_matrix(mat,0,0);
	assert(rvbool==false);
	rvbool = within_core_matrix(mat,4,0);
	assert(rvbool==false);
	rvbool = within_core_matrix(mat,1,1);
	assert(rvbool==true);
	delete_matrix(&mat);

	/* Create blinker to test conways game of life */
	mat = new_matrix(5,5);
	new_mat = new_matrix(5,5);
	set_elem(mat,2,2,1.0);

	delete_matrix(&mat);
	delete_matrix(&new_mat);		

	printf("Testing new_matrix_checkerboard_elem\n");
	checker = new_matrix_checkerboard_elem(0,4,0,0,0,0);	
	assert(checker==NULL);
	checker = new_matrix_checkerboard_elem(3,0,0,0,0,0);	
	assert(checker==NULL);
	checker = new_matrix_checkerboard_elem(3,4,-1,0,0,0);	
	assert(checker==NULL);
	checker = new_matrix_checkerboard_elem(3,4,0,-1,0,0);	
	assert(checker==NULL);
	checker = new_matrix_checkerboard_elem(3,4,0,0,-1,0);	
	assert(checker==NULL);
	checker = new_matrix_checkerboard_elem(3,4,0,0,0,-1);	
	assert(checker==NULL);
	checker = new_matrix_checkerboard_elem(3,4,0,0,0,0);	
	assert(checker!=NULL);	
	assert(get_rows(checker)==3);
	assert(get_cols(checker)==4);	
	assert(get_rows_core_matrix(checker)==3);
	assert(get_cols_core_matrix(checker)==4);
	
	printf("Testing get_cols_left_ghost_matrix & get_cols_right_ghost_matrix\n");
	rv = get_cols_left_ghost_matrix(NULL);
	assert(rv==-1);
	rv = get_cols_right_ghost_matrix(NULL);
	assert(rv==-1);
	rv = get_cols_left_ghost_matrix(checker);
	assert(rv==0);
	rv = get_cols_right_ghost_matrix(checker);
	assert(rv==0);
	delete_matrix(&checker);
	checker = new_matrix_checkerboard_elem(3,4,0,0,1,0);
	assert(checker!=NULL);
	assert(get_rows(checker)==3);
	assert(get_cols(checker)==5);
	assert(get_rows_core_matrix(checker)==3);
	assert(get_cols_core_matrix(checker)==4);
	assert(get_cols_left_ghost_matrix(checker)==0);
	assert(get_cols_right_ghost_matrix(checker)==1);
	delete_matrix(&checker);
	checker = new_matrix_checkerboard_elem(4,2,1,0,2,2);
	assert(checker!=NULL);
	print_matrix(checker);
	assert(get_rows(checker)==5);
	assert(get_cols(checker)==6);
	assert(get_rows_core_matrix(checker)==4);
	assert(get_cols_core_matrix(checker)==2);
	assert(get_cols_left_ghost_matrix(checker)==2);
	assert(get_cols_right_ghost_matrix(checker)==2);
	delete_matrix(&checker);
	return 0;
}
