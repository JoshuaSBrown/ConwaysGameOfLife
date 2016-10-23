/**************************************/
/*          Midterm 3
 *          Joshua Brown
 *          Scientific Computing
 *          midterm_matrix.h function
 *          Conway's Game of Life
 */
/**************************************/

#include <stdbool.h>

#ifndef _MIDTERM_MATRIX_H_
#define _MIDTERM_MATRIX_H_

typedef struct _matrix * matrix;
typedef struct _matrix const * const_matrix;

/* Create a new matrix of intergers */
matrix new_matrix(int rows, int cols);

/* Create a new matrix with ghost rows. Basically
 * creates a matrix with optionally ghost rows
 * above and or below. */
matrix new_matrix_row_panel(int core_rows      , int cols,
						    int ghost_row_above, int ghost_row_below);

/* Creates a new matrix with optional ghost rows
 * below, above, left and right of the core matrix
 */
matrix new_matrix_checkerboard_elem(int core_rows, int core_cols,
                                    int ghost_row_above, int ghost_row_below,
                                    int ghost_col_right, int ghost_col_left);

/* free memory used when the matrix was created */
int delete_matrix(matrix * mat);

/* matrix multiplication is executed such that the
 * output is stored in c. However, the elements in 
 * c are not overwritted but the elements are added
 * to such that you get:
 * c = c+a*b */
int matrix_multiply_and_add_to_c(matrix c, matrix a, matrix b);

/* matrix multiplication is executed such that the
 * output is stored in c. i.e.
 * c = a*b */
int matrix_multiply(matrix c, matrix a, matrix b);

/* Print the contents of a matrix                */
int print_matrix(const_matrix mat);
	
/* Print the contents of a matrix in compressed
 * form. Used to convert floating point numbers 
 * to an interger and then prints intergers with
 * place holders corresponding to two significant
 * digits */
int print_matrix_compressed(const_matrix mat);

/* Gets the start of the matrix indicates the row
 * the core matrix starts at. Any row less than 
 * than the value returned exists within the ghost
 * row on the top side                          */
int get_row_start(const_matrix mat);

/* Gets the end row of the core matrix any row
 * of equal or greater value is part of the bottom
 * ghost row.                                   */
int get_row_end(const_matrix mat);

/* Gets the col at the start of the core matrix
 * and value less than the column returned is part
 * of the left ghost column.                     */ 
int get_col_start(const_matrix mat);

/* Gets the col at the end of the core matrix any 
 * column with a value greater or equal to the 
 * value returned is part of the ghost column on
 * the right side of the core matrix.            */
int get_col_end(const_matrix mat);

/* Gets the total rows of the matrix including 
 * rows in both top and bottom ghost matrices    */
int get_rows(const_matrix mat);

/* Gets the total columns of the matrix including
 * columns in both the left and right ghost 
 * columns.                                      */ 
int get_cols(const_matrix mat);

/* Gets the total number of elements in the matrix
 * including ghost columns and rows.             */
int get_total_elems(const_matrix mat);

/* Determines the total number of elements in the 
 * core matrix. Excludes elements located in the
 * ghost rows and columns.                       */
int get_total_elems_core_matrix(const_matrix mat);

/* Gets the total number of rows in the core
 * matrix                                        */
int get_rows_core_matrix(const_matrix mat);

/* Gets the total number of cols in the core  
 * matrix                                       */
int get_cols_core_matrix(const_matrix mat);

int get_cols_left_ghost_matrix(const_matrix mat);

int get_cols_right_ghost_matrix(const_matrix mat);

int get_rows_above_ghost_matrix(const_matrix mat);

int get_rows_below_ghost_matrix(const_matrix mat);
/* Copies the elements between and including 
 * start_row to end_row
 * start_col to end_col 
 * to a separate matrix the values are placed in 
 * the new matrix starting at row and col (0,0) */
int copy_matrix_elems(int start_row_elem   ,
                      int end_row_elem     ,
                      int start_col_elem   ,
                      int end_col_elem     ,
                      const_matrix mat_orig,
                      matrix * mat_copy    );

/* Grab the element located in the global matrix
 * or full matrix this includes both the core 
 * matrix and the ghost matrices that pad the 
 * sides if any exists. The row and col refer to
 * the position within this full matrix.         */
float get_elem(const_matrix mat, int row, int col);

/* Set the element within the full or global 
 * matrix. This includes both the core matrix and
 * the ghost matrices that may bad the core 
 * matrix                                        */
int set_elem(matrix mat, int row, int col, float val);

/* Set all of the elements in the matrix mat equal
 * to val */
int set_all_elems(matrix mat, float val);

/* Set the value of an element in the core matrix
 * row and col are the row and col of the matrix.
 * they can be used as if the matrix had know ghost
 * rows or columns. For instance if the matrix does
 * have ghost rows or columns element row=1 and 
 * col=1 using this function will always refer to the
 * first element of the core matrix */
int set_elem_core_matrix(matrix mat,int row, int col, float val);

/* Get an element from the core matrix this
 * excludes any elements in the ghost. Here row
 * and column refer to the rows and cols in the
 * core matrix not the global matrix.            */
float get_elem_core_matrix(const_matrix mat, int row, int col);

/* Find the total value of all of the elements in
 * the core matrix and return the value. Excludes 
 * elements in the ghost rows or columns.        */
float sum_all_core_matrix_elems(const_matrix mat);

/* Allows you to change an element in the top 
 * ghost row matrix. Here row and col refer to  
 * the rows and columns within the ghost matrix 
 * they do not refer to the row and column of the 
 * global matrix                                 */
int set_elem_top_row_ghost_matrix(matrix mat, int row, int col, float val);

/* Grab the element in the top ghost matrix      
 * specified by row and col. Where row and col are
 * not refering to the rows and columns in the 
 * global matrix.                                */
float get_elem_top_row_ghost_matrix(const_matrix mat, int row, int col);

int set_elem_bottom_row_ghost_matrix(matrix mat, int row, int col, float val);

float get_elem_bottom_row_ghost_matrix(const_matrix mat, int row, int col);

int get_ghost_cols_left(const_matrix mat);

int get_ghost_cols_right(const_matrix mat);

int get_ghost_rows_above(const_matrix mat);

int get_ghost_rows_below(const_matrix mat);

/* Send a matrix using a derived MPI datatype    */
int send_matrix(matrix mat);

/* Receive a matrix using a derived MPI datatype */
int receive_matrix(matrix mat);

/* Send and receive ghost columns through a 
 * derived MPI datatype. The total number of 
 * columns sent and received will depend on the
 * value of ghost_count. 
 */
bool send_recv_ghost_cols(matrix * new_mat);

/* Send and receive ghost rows through a derived
 * MPI datatype. The total number of rows sent
 * and recieved will depend on the value of 
 * ghost_count.                                  */
bool send_recv_ghost_rows(matrix * new_mat);

/* Determine if a cell is alive depending on the 
 * total number of bugs (sum) from the surrounding
 * cells.                                        */
float live_cell(float sum, float cell);

static inline void conways_condition(const_matrix orig_mat, 
                                     matrix * new_mat, 
                                     int i, int j,
                                     int ghost);

/* Checks to see if the elements i,j is within 
 * the core matrix if it is returns true if not
 * it returns false                              */
bool within_core_matrix(const_matrix mat, int i, int j);

bool within_core_matrix_with_ghost(const_matrix mat, 
                                   int i, int j,
                                   int ghost);

/* Checks to see if the element i,j is within the
 * matrix mat. If it is it returns true otherwise
 * returns false.                                */
bool within_matrix(const_matrix mat, int i, int j);

/* Simply generates an ascii data file for matlab
 * does not work in parallel                     */
int generate_data_file(char * filename, matrix mat);

/* This function can run conways game for both 
 * checkerboard and block row decompositions. It
 * is also able to take advantage of the ghost_rows
 * and ghost_columns if they exceed a single layer
 * this is done by doing extra calculations in
 * place of sending messages. 
 */
bool conways_game(matrix * orig_mat, 
                  matrix * temp, 
                  matrix * new_mat);

/* This function is used to test. It simply add all
 * of the elements surrouding i,j and returns the 
 * sum                                           */
float add_all_surrounding_elems(const_matrix orig_mat, 
                                int i, int j);

/* Designed to run conway's game of life on the
 * border of the matrix. The thickness of the 
 * border is determined by ghost_count. In this
 * case conway's game of life is calculated for
 * the values stored in temp1 the result is put
 * in temp2. Returns a value of 0 if it exectued
 * correctly, if not returns a -1.               */
int update_border_conway(matrix temp1, 
                         matrix * temp2,
                         int k);

/* Used to update border of matrix executed 
 * conways on mat1 and stores the result in mat2 */
int set_elems_conway_rectangle(int col_start, 
                               int col_end,
                               int row_start, 
                               int row_end,
                               const_matrix mat1, 
                               matrix * mat2,
                               int ghost);

/* Does the same thing as the above function with
 * the exception that when k is equal to the
 * (ghost_count-1) the final value is stored in 
 * new. This function is used when calculating
 * the values within the margin.                 */
int set_elems_conway_rectangle_rotation(
                               int col_start,
                               int col_end,
                               int row_start,
                               int row_end,
                               const_matrix mat1,
                               matrix *mat2,
						       matrix *new_mat, 
                               int k);
	                                   
int update_conways_right(matrix orig,
						  matrix temp, 
						  matrix new_mat,
						  int k);

/*int update_conways_right_odd(matrix orig,
                             matrix temp,
                             matrix new_mat,
                             int k);
*/
int update_conways_left(matrix orig,
						 matrix temp,
						 matrix new_mat,
						 int k);
/*
int update_conways_left_odd(matrix orig,
                            matrix temp,
                            matrix new_mat,
                            int k);
*/
#endif
