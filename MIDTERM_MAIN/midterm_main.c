/**************************************/
/*          Midterm 3
 *          Joshua Brown
 *          Scientific Computing
 *          midterm_main.c function
 *          Conway's Game of Life
 */
/**************************************/

#include <stdio.h>
#include <stdlib.h>


#include "../MIDTERM_MATRIX/midterm_matrix.h"
#include "../MIDTERM_IO/midterm_io.h"
#include "../MIDTERM_CONSTANTS/midterm_constants.h"
#include "../MIDTERM_ERROR/midterm_error.h"

#ifdef _MPI_H_
#include "mpi.h"
#endif 
int main(int argc, char **argv){

	#ifdef _MPI_H_
	MPI_Init(0,0);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&num_proc);
	#else
	my_rank = 0;
	num_proc = 1;
	#endif
	
	char timing_file[50];
	/* Creating and open timing file if rank 0 */
	if(my_rank==0){
		sprintf(timing_file,"timing.txt");
	}
	printf("Timer resolution %g ",MPI_Wtick());
	stat_sum = 0.0;
	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();

	float local_bug_sum;
	int counter;
	int pgm_increment_counter;
	int i;
	int j;
	int k;
	float val;
	double total_time_core = 0.0;
	unsigned long int bytes_core;
	unsigned long int bytes_border;
	unsigned long int bytes_message;
	int ghost_left;
	int ghost_right;
	int ghost_above;
	int ghost_below;

	create_filename_out(my_rank);
	switch_stdout(filename_out);	
	/* Parse input arguments */
	options(argc, argv);

	/* Initialize the following global variables:
     *    o grid_dim_col
     *    o grid_dim_row
     *    o my_grid_col
     *    o my_grid_row
	 */
	if(block_on){
		initialize_grid_row_block();
	}else{
		initialize_grid_checker();
	}
	matrix mat_a = NULL;
	matrix mat_b = NULL;
	matrix extra = NULL;
	matrix * temp1 = NULL;
	matrix * temp2 = NULL;
	matrix * temp3 = NULL;
	pgm_increment_counter = 0;

	if(init_test==0){
		read_pgm(&mat_a, &mat_b);
	}else{
		int width;
		int height;
		width  = grid_dim_col*12;
		height = grid_dim_row*12;
		sub_matrix_dim_columns = width  / grid_dim_col;
		sub_matrix_dim_rows    = height / grid_dim_row;
		mat_a = define_matrix_blockchecker();
		mat_b = define_matrix_blockchecker();
		
		/* We will begin by setting all the elements
         * to 1's                                  */
		set_all_elems(mat_a,1.0);
		set_all_elems(mat_b,1.0);
	}

	extra = define_matrix_blockchecker();
	
	temp1 = &mat_a;
	temp3 = &extra; 
	/* Initialize the number of bugs in the system */
	global_bug_sum = 0;	
	local_bug_sum  = 0;
	/***********************************************************************/
	/* CREATE FIRST .pgm FILE                                              */
	/***********************************************************************/
	if(init_matrix_increment_pgm!=NULL){
		write_pgm(*temp1, pgm_increment_counter); 		
		pgm_increment_counter++;	
	}
	/***********************************************************************/
	/* END OF CREATING FIRST .pgm FILE                                     */
	/***********************************************************************/

	/* Determine how many bytes will be sent via message */
	bytes_message = 0;
	bytes_border  = 0;	
	bytes_core    = 0;

        ghost_left  = get_cols_left_ghost_matrix(mat_a);
        ghost_right = get_cols_right_ghost_matrix(mat_a);
        ghost_above = get_rows_above_ghost_matrix(mat_a);
        ghost_below = get_rows_below_ghost_matrix(mat_a);
	int left  = 0;
	int above = 0;
	int below = 0;
	int right = 0;
	if(num_proc==1){
		bytes_message = 0;
		bytes_border  = 0;
		bytes_core    = get_total_elems_core_matrix(mat_a);
		
	}else{
		if(block_on){

			if(my_rank==0 || my_rank==(num_proc-1)){
				bytes_message += (ghost_count)*get_cols_core_matrix(mat_a);
				bytes_border  += get_cols_core_matrix(mat_a)*ghost_count;
				bytes_core    += get_cols_core_matrix(mat_a)*(get_rows_core_matrix(mat_a)-ghost_count);	
			}else{
				bytes_message += 2*(ghost_count)*get_cols_core_matrix(mat_a);
				bytes_border  += 2*get_cols_core_matrix(mat_a)*ghost_count;	
				bytes_core    += get_cols_core_matrix(mat_a)*(get_rows_core_matrix(mat_a)-2*ghost_count);	
			}

			for(counter=1;counter<ghost_count;counter++){
				if(my_rank==0 || my_rank==(num_proc-1)){
					bytes_border += get_cols_core_matrix(mat_a)*(ghost_count+2*counter);
					bytes_core   += get_cols_core_matrix(mat_a)*(get_rows_core_matrix(mat_a)-ghost_count+counter);	
				}else{
					bytes_border += 2*get_cols_core_matrix(mat_a)*(ghost_count+2*counter);
					bytes_core   += get_cols_core_matrix(mat_a)*(get_rows_core_matrix(mat_a)-2*ghost_count+2*counter);	
				}
			}
		}else{

			bytes_border = (get_total_elems_core_matrix(mat_a))-
				(get_rows_core_matrix(mat_a)-ghost_above-ghost_below)*
				(get_cols_core_matrix(mat_a)-ghost_left-ghost_right);
			bytes_message = bytes_border;

			bytes_core = get_total_elems_core_matrix(mat_a)-bytes_border;

			for(counter=1;counter<ghost_count;counter++){

				if(ghost_above!=0){
					above = counter;
				}
				if(ghost_below!=0){
					below = counter;
				}
				if(ghost_left!=0){
					left = counter;
				}
				if(ghost_right!=0){
					right = counter;
				}

				bytes_border += ((get_rows_core_matrix(mat_a)+above+below)*
						(get_cols_core_matrix(mat_a)+left+right))-
					(get_rows_core_matrix(mat_a)-ghost_above-ghost_below-above-below)*
					(get_cols_core_matrix(mat_a)-ghost_left-ghost_right-left-right);

				bytes_core   += (get_rows_core_matrix(mat_a)-ghost_above-ghost_below+above+below)*
					(get_cols_core_matrix(mat_a)-ghost_left-ghost_right+left+right);
			}

		}
	}
	/***********************************************************************/
	/* INITIALIZE GHOST ROWS AND COLUMNS BY COMMUNICATING WITH NEIGHBORS   */
	/***********************************************************************/
	if(!(my_rank==(num_proc-1) && my_rank==0) ){

		/******************************************************************/
		/* FOR CHECKER ONLY                                               */
		/******************************************************************/

		if(block_on==false){
			if(verbose){
				printf("send_recv_ghost_cols at counter %d\n",0);
			}
			send_recv_ghost_cols(temp1);

			if(send_recv_sync==false){
				/* Asynchronous comminication is turned on */
				if(verbose){
					printf("MPI_Wait for cols at counter request3"
							" %d\n",0);
				}
				#ifdef _MPI_H_
				if(my_grid_col==0){
					MPI_Wait(&send_request3,&status);
				}else if (my_grid_col==grid_dim_col-1){
					MPI_Wait(&recv_request3,&status);
				}else{	
					MPI_Wait(&recv_request3,&status);
					MPI_Wait(&send_request3,&status);
				}
				#endif
				if(verbose){
					printf("MPI_Wait for cols at counter request4"
							" %d\n",0);
				}
				#ifdef _MPI_H_
				if(my_grid_col==0){
					MPI_Wait(&recv_request4,&status);
				}else if (my_grid_col==grid_dim_col-1){
					MPI_Wait(&send_request4,&status);
				}else{	
					MPI_Wait(&recv_request4,&status);
					MPI_Wait(&send_request4,&status);
				}
				#endif

			}
		}
		/******************************************************************/
		/* END OF CHECKER ONLY                                            */
		/******************************************************************/

		/******************************************************************/
		/* ROW BLOCK AND CHECKER                                          */
		/******************************************************************/

		send_recv_ghost_rows(temp1);
		if(send_recv_sync==false){
			#ifdef _MPI_H_
			if(my_grid_row==0){
				MPI_Wait(&send_request1,&status);
			}else if (my_grid_row==grid_dim_row-1){
				MPI_Wait(&recv_request1,&status);
			}else{	
				MPI_Wait(&recv_request1,&status);
				MPI_Wait(&send_request1,&status);
			}
			if(my_grid_row==0){
				MPI_Wait(&recv_request2,&status);
			}else if (my_grid_row==grid_dim_row-1){
				MPI_Wait(&send_request2,&status);
			}else{	
				MPI_Wait(&recv_request2,&status);
				MPI_Wait(&send_request2,&status);
			}
			#endif
		}	
		/******************************************************************/
		/* END OF ROW BLOCK AND CHECKER                                   */
		/******************************************************************/
	}


	finish = MPI_Wtime();
	elapsed = (finish-start);
	calculate_statistics(-1);
	if(my_rank==0){
		printf("Initial Startup %g\n",max_t);
		printf("Message size %ld per %d iterations\n"
                       ,bytes_message*sizeof(float),ghost_count);
		printf("Number of floats operated on within core"
                       " of my world %ld per %d iterations\n"
                       ,bytes_core*sizeof(float),ghost_count);
		printf("Number of floats operated on within border"
                       " of my world %ld per %d iterations\n"
                       ,bytes_border*sizeof(float),ghost_count);
		printf("Total number of floats operated on %ld per"
                       " %d interations\n"
                       ,(bytes_border+bytes_core)*sizeof(float)
                       ,ghost_count);
	}
	stat_sum = 0.0;
	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();
		
	/* Begin Conway's Game of Life */
	for( counter=0; counter<init_iter; counter++){

		/* This means we want to dump output need to check if increment is
         * of interest to user */
		if(counter%ghost_count==0){
		/**************************************************************/
		/* COUNT AND PRINT BUG SUM                                    */
		/**************************************************************/
		if(counter==0){
			local_bug_sum = sum_all_core_matrix_elems(*temp1);
			if(verbose){
				printf("Calling MPI_Reduce at counter %d\n",counter);
			}
			#ifdef _MPI_H_
			MPI_Reduce(&local_bug_sum,&global_bug_sum,1,
					MPI_FLOAT,MPI_SUM, 0,MPI_COMM_WORLD);
			#endif
			if(my_rank==0){
				printf("\nIteration %6d Global sum %6d local sum %6d\n",
						counter, (int)global_bug_sum,(int)local_bug_sum);
				fflush(NULL);
			}else{
				printf("\nIteration %6d Local sum %6d\n",
                       counter, (int)local_bug_sum);
				fflush(NULL);
			}
		}else if(init_print>0){
			if(counter%init_print==0){
				local_bug_sum = sum_all_core_matrix_elems(*temp1);
				#ifdef _MPI_H_
				MPI_Reduce(&local_bug_sum,&global_bug_sum,1,
						MPI_FLOAT,MPI_SUM, 0,MPI_COMM_WORLD);
				#endif
				if(my_rank==0){
					printf("Iteration %6d Global sum %6d Local sum %6d\n",
							counter,(int) global_bug_sum, (int) local_bug_sum);
					fflush(NULL);
				}else{
					printf("Iteration %6d Local sum %6d\n",
                            counter, (int) local_bug_sum);
				}
			}
		}
		/****************************************************************/
		/* END OF PRINT BUG SUM                                         */
		/****************************************************************/

		/****************************************************************/
		/* CONWAY'S GAME OF LIFE BEGIN                                  */
		/****************************************************************/
		conways_game(temp1,&mat_b, temp3);
		/****************************************************************/
		/* CONWAY'S GAME OF LIFE END                                    */
		/****************************************************************/
		/* Alternative which matrix is considered new 
		 * After conways temp3 will always be the newest */
		temp2 = temp1;
		temp1 = temp3;
		temp3 = temp2;

		/* This means we want to dump output need to check if increment is
		 * of interest to user */

		/****************************************************************/
        /* FILE PRINTER                                                 */
		/****************************************************************/
		if(init_matrix_increment_pgm!=NULL){
			if(flag_increment_pgm==2){
				/* Interested in dumping .pgm file if appears within range */
				if(counter>=get_elem(init_matrix_increment_pgm,0,0) && 
						counter<=get_elem(init_matrix_increment_pgm,1,0)){
					write_pgm(*temp1,pgm_increment_counter); 		
					pgm_increment_counter++;	
				}

			}else if(flag_increment_pgm==1){
				/* Interested in dumping .pgm file if appears in list */			
				for(k=0;k<get_rows(init_matrix_increment_pgm);k++){
					if(counter==get_elem(init_matrix_increment_pgm,k,0)){
						write_pgm(*temp1,pgm_increment_counter); 			
						pgm_increment_counter++;	
						break;
					}
				}
			}else if(flag_increment_pgm==0){
				/* Only interested in dumping a single .pgm file */
				if(counter==get_elem(init_matrix_increment_pgm,0,0)){
					write_pgm(*temp1,pgm_increment_counter); 			
					pgm_increment_counter++;	
				}
			}
		} 
		/****************************************************************/
		/* END OF FILE PRINTER                                          */
		/****************************************************************/
		}

		finish = MPI_Wtime();
		elapsed = (finish-start);
		calculate_statistics(counter);
		total_time_core+=elapsed;
		if(my_rank==0){
			generate_performance_file(timing_file,counter,bytes_message);
		}
               	MPI_Barrier(MPI_COMM_WORLD);
		start = MPI_Wtime();

	}/* Counter loop */	

	finish = MPI_Wtime();
	elapsed = (finish-start);
	calculate_statistics(1);
	printf("Program Core    %g Total iterations %d"
               " Ghost count     %d\n"
               ,total_time_core,init_iter,ghost_count);
	stat_sum = 0.0;
	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();

	local_bug_sum = sum_all_core_matrix_elems(*temp1);
	#ifdef _MPI_H_
	MPI_Reduce(&local_bug_sum,&global_bug_sum,1,MPI_FLOAT,
               MPI_SUM, 0,MPI_COMM_WORLD);
	#endif
	if(my_rank==0){
		printf("Iteration %6d Global sum %6d\n",counter, (int)global_bug_sum);
	}
	if(init_matrix_increment_pgm!=NULL){
		delete_matrix(&init_matrix_increment_pgm);
	}

	/* Delete Matrices */
	delete_matrix(&mat_a);
	delete_matrix(&mat_b);
	delete_matrix(&extra);
	revert_stdout();
	
	finish = MPI_Wtime();
	elapsed = (finish-start);
	calculate_statistics(1);
        printf("Program End     %g\n",max_t);
	#ifdef _MPI_H_
	MPI_Finalize();
	#endif	

	return 0;
}
