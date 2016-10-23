/**************************************/
/*          Midterm 3
 *          Joshua Brown
 *          Scientific Computing
 *          midterm_constants.h function
 *          Conway's Game of Life
 */
/**************************************/
#include <stdbool.h>
#include "../MIDTERM_ERROR/midterm_error.h"
#include "../MIDTERM_MATRIX/midterm_matrix.h"

#ifndef _MIDTERM_CONSTANTS_
#define _MIDTERM_CONSTANTS_

#ifdef _MPI_H_
#include "mpi.h"
#endif
/* Processor rank                           */
int my_rank;

/* Number of Processors                     */
int num_proc;

/* What is the grid location of my_rank     */
int my_grid_row;
int my_grid_col;

/* What are the dimension of the grid in terms
 * of all the processors                    */
int grid_dim_row;
int grid_dim_col;

/* rank of processor in a row of procs      */
int my_row_rank;

/* Number of processors in the row          */
int num_proc_in_row;

/* Communicator for the row                 */
#ifdef _MPI_H_
MPI_Comm row_comm;

/* MPI status used when receiving messages  */
MPI_Status status;

/* Used with both checkerboard and panel
 * decomposition.                           */
MPI_Request send_request1, recv_request1;
MPI_Request send_request2, recv_request2;

/* To be used exclusively with checkerboard */
MPI_Request send_request3, recv_request3;
MPI_Request send_request4, recv_request4;
#endif

/* Size of the sub_matrix dimension         */
int sub_matrix_dim_columns;
int sub_matrix_dim_rows;

/* Source rank where we are receiving a 
 * message from.                            */
int source;

/* Destination rank where we are sending a 
 * message too                              */
int dest;

/* Destination rank where we are sending the
 * message                                  */
int message;

/* Global sum of bugs reduced from all  
 * processors.                              */
float global_bug_sum;

/* Options from user                        */
bool   verbose;
bool   block_on;
bool   send_recv_sync;
int    init_rows;
int    init_cols;
int    init_iter;
int    init_print;
bool   init_data;
char * init_print_pgm;
int    init_test;
char * init_file;
int    ghost_count;

/* This matrix stores the values that are 
 * meant to be printed out                  */
matrix init_matrix_increment_pgm;
int flag_increment_pgm;

/* Performance variables                    */
double tick;
double start;
double finish;
double elapsed;
double max_t;
double stat_sum;
double running_mean;
double standard_dev;
double margin_err;
double t_val;

/* Determine constants:
 *   my_grid_row
 *   my_grid_col
 *   grid_dim_row
 *   grid_dim_col
 * For checker decomposition                */
int initialize_grid_checker(void);

/* Determine constants:
 *   my_grid_row
 *   my_grid_col
 *   grid_dim_row
 *   grid_dim_col
 * For panel decomposition                  */
int initialize_grid_row_block(void);

/* Given a grid col and row will determine  
 * what rank/index should be associated with 
 * it.                                      */
int grid_index(int friend_grid_col, 
               int friend_grid_row);
#endif
