/**************************************/
/*          Midterm 2
 *          Joshua Brown
 *          Scientific Computing
 *          midterm_constants.c function
 *          Conway's Game of Life
 */
/**************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "midterm_constants.h"

int initialize_grid_checker(void){

	float sqrt = pow((float)num_proc,1/2.0);

	if(fmod(sqrt,1.0)==0.0){
		grid_dim_col = (int) sqrt;
		grid_dim_row= (int) sqrt;

		my_grid_col= my_rank%grid_dim_row;
		my_grid_row= my_rank/grid_dim_row;

		printf("In function my_rank %d my_grid_col %d my_grid_row %d\n",my_rank,my_grid_col,my_grid_row);
	}else{
		fprintf(stderr,"ERROR num of processors does not divide evenly up into"
                       " columns and rows int initialize_grid_checker\n");
		exit(1);
	}

	return 0;
}

int initialize_grid_row_block(void){

	my_grid_col = 0;
	my_grid_row = my_rank; 
	grid_dim_col = 1;
	grid_dim_row = num_proc;
	
	return 0;
}

int grid_index(int friend_grid_row, int friend_grid_col){
	
	return friend_grid_row*(grid_dim_col)+friend_grid_col;	
}
