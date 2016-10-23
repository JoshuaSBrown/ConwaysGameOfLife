/**************************************/
/*          Midterm 3 
 *          Joshua Brown
 *          Scientific Computing
 *          midterm_matrix.c function
 *          Conway's Game of Life
 */
/**************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "midterm_matrix.h"
#include "../MIDTERM_ERROR/midterm_error.h"
#include "../MIDTERM_CONSTANTS/midterm_constants.h"

#ifdef _MPI_H_
#include "mpi.h"
#endif

#define ELEM(mat, row, col) \
  ((mat)->data[(row) * (mat)->cols + (col)])

/* Here rows and cols describe the height and width
 * of the entire matrix. However, it is possible to
 * also have ghost rows and columns all around the 
 * core matrix. Thus we have defined for parmeters
 * to define where the core matrix begins and end 
 * in terms of the cols and the rows. */
struct _matrix {
	int   rows;
	int   cols;
	int   col_start;
	int   col_end;
	int   row_start;
	int   row_end;
	float data[0];
};

/* Creates a 2d matrix from scratch */
matrix new_matrix(int rows, int cols){

	#ifdef _ERROR_CHECKING_ON_
	if(rows<1){
		fprintf(stderr,"ERROR rows is less than 1 in new_matrix\n");
		return_error_val();
		return NULL;
	}
	if(cols<1){
		fprintf(stderr,"ERROR cols is less than 1 in new_matrix\n");
		return_error_val();
		return NULL;
	}
	#endif

	matrix mat = (matrix) malloc(sizeof(struct _matrix)+sizeof(float)*(rows+1)*cols);

	#ifdef _ERROR_CHECKING_ON_
	if(mat==NULL){
		fprintf(stderr,"ERROR mat is NULL in new_matrix\n");
		return_error_val();
		return NULL;
	}
	#endif
	
	mat->rows=rows;
	mat->cols=cols;
	mat->col_start = 0;
	mat->col_end   = cols;
	mat->row_start = 0;
	mat->row_end   = rows;	
	/* Initialize matrix so each of the elements 
     * is 0 */
	int i;
	int j;
	for(i=0;i<rows;i++){
		for(j=0;j<cols;j++){
			ELEM(mat,i,j)=0;
		}
	}

	return mat;
}

/* Creates a 2d matrix from scratch, with ghost extensions 
 * the extenions can extend above and below the matrix defined
 * by rows and cols, the ghost_rows above and below the core
 * matrix must be of the same length if they are not set to 0.
 */
matrix new_matrix_row_panel(int core_rows	   , int cols, 
							int ghost_row_above, int ghost_row_below ){

	#ifdef _ERROR_CHECKING_ON_
	if(core_rows<1){
		fprintf(stderr,"ERROR core_rows is less than 1 in "
                       "new_matrix_row_panel\n");
		return_error_val();
		return NULL;
	}
	if(cols<1){
		fprintf(stderr,"ERROR cols is less than 1 in "
                       "new_matrix_row_panel\n");
		return_error_val();
		return NULL;
	}
	if(ghost_row_above<0){
		fprintf(stderr,"ERROR ghost_row_above is less than 0 "
                       "in new_matrix_row_panel\n");
		return_error_val();
		return NULL;
	}
	if(ghost_row_below<0){
		fprintf(stderr,"ERROR ghost_row_below is less than 0 "
                       "in new_matrix_row_panel\n");
		return_error_val();
		return NULL;
	}
	#endif

	int total_rows = core_rows+ghost_row_above+ghost_row_below;

	matrix mat = (matrix) malloc(sizeof(struct _matrix)+sizeof(float)*(total_rows)*cols);

	#ifdef _ERROR_CHECKING_ON_
	if(mat==NULL){
		fprintf(stderr,"ERROR mat is NULL in new_matrix\n");
		return_error_val();
		return NULL;
	}
	#endif
	
	mat->rows=total_rows;
	mat->cols=cols;
	mat->row_start = ghost_row_above;
	mat->row_end   = core_rows+ghost_row_above;
	mat->col_start = 0;
	mat->col_end   = cols;
	/* Initialize matrix so each of the elements 
     * is 0 */
	int i;
	int j;
	for(i=0;i<mat->rows;i++){
		for(j=0;j<cols;j++){
			ELEM(mat,i,j)=0;
		}
	}

	return mat;
}

/* Creates a 2d matrix from scratch, with ghost extensions 
 * the extenions can extend above and below the matrix defined
 * by rows and cols, they can also extend to the left and the
 * right. The ghost_rows above and below the core matrix must
 * be of the same length unless they are of length 0. 
 * The same goes for the ghost columns */
matrix new_matrix_checkerboard_elem(int core_rows	   , int core_cols, 
							        int ghost_row_above, int ghost_row_below, 
                                    int ghost_col_right, int ghost_col_left){

	#ifdef _ERROR_CHECKING_ON_
	if(core_rows<1){
		fprintf(stderr,"ERROR core_rows is less than 1 in " 
                       "new_matrix_checkerboard_elem\n");
		return_error_val();
		return NULL;
	}
	if(core_cols<1){
		fprintf(stderr,"ERROR cols is less than 1 in "
                       "new_matrix_checkerboard_elem\n");
		return_error_val();
		return NULL;
	}
	if(ghost_row_above<0){
		fprintf(stderr,"error ghost_row_above is less than 0 "
                       "in new_matrix_checkerboard_elem\n");
		return_error_val();
		return NULL;
	}
	if(ghost_row_below<0){
		fprintf(stderr,"error ghost_row_below is less than 0 "
                       "in new_matrix_checkerboard_elem\n");
		return_error_val();
		return NULL;
	}
	if(ghost_col_right<0){
		fprintf(stderr,"error ghost_col_right is less than 0 "
                       "in new_matrix_checkerboard_elem\n");
		return_error_val();
		return NULL;
	}
	if(ghost_col_left<0){
		fprintf(stderr,"error ghost_col_left is less than 0 in"
                       " new_matrix_checkerboard_elem\n");
		return_error_val();
		return NULL;
	}
	#endif

	int total_rows = core_rows+ghost_row_above+ghost_row_below;
	int total_cols = core_cols+ghost_col_left +ghost_col_right;

	matrix mat = (matrix) malloc(sizeof(struct _matrix)+sizeof(float)*(total_rows)*total_cols);

	#ifdef _ERROR_CHECKING_ON_
	if(mat==NULL){
		fprintf(stderr,"ERROR mat is NULL in "
                       "new_matrix_checkerboard_elem\n");
		return_error_val();
		return NULL;
	}
	#endif
	
	mat->rows=total_rows;
	mat->cols=total_cols;
	mat->row_start = ghost_row_above;
	mat->row_end   = core_rows+ghost_row_above;
	mat->col_start = ghost_col_left;
	mat->col_end   = core_cols+ghost_col_left;
	/* Initialize matrix so each of the elements 
     * is 0 */
	int i;
	int j;
	for(i=0;i<mat->rows;i++){
		for(j=0;j<mat->cols;j++){
			ELEM(mat,i,j)=0;
		}
	}

	return mat;
}

int delete_matrix(matrix * mat){

	#ifdef _ERROR_CHECKING_ON_
	if(mat==NULL){
		fprintf(stderr,"ERROR mat is NULL in delete_matrix\n");
		return return_error_val();
	}
	#endif

	free(*mat);
	*mat = NULL;
	return 0;
}

/* Function multiplies two matrices together to get a third matrix 
 * then adds to what ever values were already stored in matrix_c:
 * c = c+a*b
 */
int matrix_multiply_and_add_to_c(matrix c, matrix a, matrix b){

	#ifdef _ERROR_CHECKING_ON_	
	if(a->cols!=b->rows){ 
		fprintf(stderr,"ERROR a matrix does not have the "
                       "number of columns as b has rows "
                       "in matrix_multiply_and_add_to_c\n");
		return return_error_val();
	}
	if(a->rows!=c->rows){
		fprintf(stderr,"ERROR c matrix does not have the "
                       "number of rows as a in "
                       "matrix_multiply_and_add_to_c\n");
		return return_error_val();
	}
	if(b->cols!=c->cols){
		fprintf(stderr,"ERROR c matrix does not have the "
                       "number of cols as b in "
                       "matrix_multiply_and_add_to_c\n");
		return return_error_val();
	}
	#endif

	int i, j, k;
	float val;	
	for(k=0;k<a->cols;k++){
		for(i=0;i<a->rows;i++){
			val = 0.0;
			for(j=0;j<a->cols;j++){
				val+=ELEM(a,i,j)*ELEM(b,j,k);
			}
			ELEM(c,i,k) += val;
		}
	}

	return 0;
}

/* Function multiplies two matrices together to get a third matrix 
 * Such that:
 * c = a*b
 */
int matrix_multiply(matrix c, matrix a, matrix b){

	#ifdef _ERROR_CHECKING_ON_	
	if(a->cols!=b->rows){ 
		fprintf(stderr,"ERROR a matrix does not have the "
                       "number of columns as b has rows "
                       "in matrix_multiply\n");
		return return_error_val();
	}
	if(a->rows!=c->rows){
		fprintf(stderr,"ERROR c matrix does not have the "
                       "number of rows as a in matrix_multiply\n");
		return return_error_val();
	}
	if(b->cols!=c->cols){
		fprintf(stderr,"ERROR c matrix does not have the "
                       "number of cols as b in "
                       "matrix_multiply\n");
		return return_error_val();
	}
	#endif

	int i, j, k;
	float val;	

	for(k=0;k<a->cols;k++){
		for(i=0;i<a->rows;i++){
			val = 0.0;
			for(j=0;j<a->cols;j++){
				val+=ELEM(a,i,j)*ELEM(b,j,k);
			}
			ELEM(c,i,k) = val;
		}
	}

	fprintf(stderr,"Matrix Multiplication \n\n");

	return 0;
}

int print_matrix(const_matrix mat){

	int i;
	int j;
	#ifdef _ERROR_CHECKING_ON_
	if(mat==NULL){
		fprintf(stderr,"ERROR mat is NULL in print_matrix\n");
		return return_error_val();
	}
	#endif
	printf("\n*****************************************\n");
	printf("Total cols %d Total rows %d\n\n",mat->cols,mat->rows);
	printf("col_start %d col_end %d\n",mat->col_start,mat->col_end);
	printf("row_start %d row_end %d\n",mat->row_start,mat->row_end);	

	if(mat->row_start>0){
		printf("Ghost Rows %d\n",mat->row_start);
		i = 0;
		for(j=0;j<mat->cols;j++){
			printf("\tCol %d\t",j+1);	
		}
		printf("\n");

		for(i=0;i<(mat->row_start);i++){
			printf("Row %d\t",i+1);
			for(j=0;j<mat->cols;j++){
				printf("%6f\t",ELEM(mat,i,j));
			}
		}
		printf("\n");
	}

	printf("Core Matrix Rows %d\n",mat->row_end-mat->row_start);
	for(j=0;j<mat->cols;j++){
		if(j<mat->col_start){
			printf("\t Ghost Col %d",j+1);
		}else if(j>=mat->col_end){
			printf("\t Ghost Col %d",j+1);
		}else{
			printf("\t Col %d\t",j+1);
		}	
	}
	printf("\n");

	for(i=(mat->row_start);i<mat->row_end;i++){
		printf("Row %d\t",i+1);
		for(j=0;j<mat->cols;j++){
			printf("%6f\t",ELEM(mat,i,j));
		}
		printf("\n");
	}
	printf("\n");

	if(mat->row_end!=mat->rows){
		printf("Ghost Rows %d\n",(mat->rows-mat->row_end));
		for(j=0;j<mat->cols;j++){
			printf("\t Col %d\t",j+1);	
		}
		printf("\n");
		for(i = mat->row_end;i<mat->rows;i++){
			printf("Row %d\t",i+1);
			for(j=0;j<mat->cols;j++){
				printf("%6f\t",ELEM(mat,i,j));
			}
			printf("\n");
		}
		printf("\n");
	}

	return 0;
}

int print_matrix_compressed(const_matrix mat){

	#ifdef _ERROR_CHECKING_ON_
	if(mat==NULL){
		fprintf(stderr,"ERROR mat is NULL in print_matrix\n");
		return return_error_val();
	}
	#endif
	int i;
	int j;
	printf("\n*****************************************\n");
	printf("Total cols %d Total rows %d\n\n",mat->cols,mat->rows);
	printf("col_start %d col_end %d\n",mat->col_start,mat->col_end);
	printf("row_start %d row_end %d\n",mat->row_start,mat->row_end);	

	if(mat->row_start>0){
		for(i=0;i<(mat->row_start);i++){
			for(j=0;j<mat->cols;j++){
				printf("%3d ",(int)ELEM(mat,i,j));
			}
			printf("\n");
		}
	}

	for(i=(mat->row_start);i<mat->row_end;i++){
		for(j=0;j<mat->cols;j++){
			printf("%3d ",(int)ELEM(mat,i,j));
		}
		printf("\n");
	}

	if(mat->row_end!=mat->rows){
		for(i = mat->row_end;i<mat->rows;i++){
			for(j=0;j<mat->cols;j++){
				printf("%3d ",(int)ELEM(mat,i,j));
			}
			printf("\n");
		}
		printf("\n");
	}

	return 0;
}

int get_row_start(const_matrix mat){
	#ifdef _ERROR_CHECKING_ON_
	if(mat==NULL){
		fprintf(stderr,"ERROR mat is NULL in get_row_start\n");
		return return_error_val();
	}
	#endif
	return mat->row_start;
}
int get_row_end(const_matrix mat){
	#ifdef _ERROR_CHECKING_ON_
	if(mat==NULL){
		fprintf(stderr,"ERROR mat is NULL in get_row_end\n");
		return return_error_val();
	}
	#endif
	return mat->row_end;
}
int get_col_start(const_matrix mat){
	#ifdef _ERROR_CHECKING_ON_
	if(mat==NULL){
		fprintf(stderr,"ERROR mat is NULL in get_col_start\n");
		return return_error_val();
	}
	#endif
	return mat->col_start;
}
int get_col_end(const_matrix mat){
	#ifdef _ERROR_CHECKING_ON_
	if(mat==NULL){
		fprintf(stderr,"ERROR mat is NULL in get_col_end\n");
		return return_error_val();
	}
	#endif
	return mat->col_end;
}

int get_rows(const_matrix mat){
	#ifdef _ERROR_CHECKING_ON_
	if(mat==NULL){
		fprintf(stderr,"ERROR mat is NULL in get_rows\n");
		return return_error_val();
	}
	#endif
	return mat->rows;
}

int get_cols(const_matrix mat){
	#ifdef _ERROR_CHECKING_ON_
	if(mat==NULL){
		fprintf(stderr,"ERROR mat is NULL in get_cols\n");
		return return_error_val();
	}
	#endif
	return mat->cols;
}

int get_total_elems(const_matrix mat){
	#ifdef _ERROR_CHECKING_ON_
	if(mat==NULL){
		fprintf(stderr,"ERROR mat is NULL in get_total_elems\n");
		return return_error_val();
	}
	#endif
	return mat->cols*mat->rows;
}

int get_total_elems_core_matrix(const_matrix mat){
	#ifdef _ERROR_CHECKING_ON_
	if(mat==NULL){
		fprintf(stderr,"ERROR mat is NULL in get_total_elems\n");
		return return_error_val();
	}
	#endif
	return get_rows_core_matrix(mat)*get_cols_core_matrix(mat);
}

int get_rows_core_matrix(const_matrix mat){
	#ifdef _ERROR_CHECKING_ON_
	if(mat==NULL){
		fprintf(stderr,"ERROR mat is NULL in get_rows_core_matrix\n");
		return return_error_val();
	}
	#endif
	return (mat->row_end-mat->row_start);
}

int get_cols_core_matrix(const_matrix mat){
	#ifdef _ERROR_CHECKING_ON_
	if(mat==NULL){
		fprintf(stderr,"ERROR mat is NULL in get_cols_core_matrix\n");
		return return_error_val();
	}
	#endif
	return (mat->col_end-mat->col_start);
}

int get_cols_left_ghost_matrix(const_matrix mat){
	#ifdef _ERROR_CHECKING_ON_
	if(mat==NULL){
		fprintf(stderr,"ERROR mat is NULL in get_cols_left_ghost_matrix\n");
		return return_error_val();
	}
	#endif
	return (mat->col_start);
}

int get_cols_right_ghost_matrix(const_matrix mat){
	#ifdef _ERROR_CHECKING_ON_
	if(mat==NULL){
		fprintf(stderr,"ERROR mat is NULL in get_cols_right_ghost_matrix\n");
		return return_error_val();
	}
	#endif
	return (mat->cols-mat->col_end);
}

int get_rows_above_ghost_matrix(const_matrix mat){
	#ifdef _ERROR_CHECKING_ON_
	if(mat==NULL){
		fprintf(stderr,"ERROR mat is NULL in get_rows_above_ghost_matrix\n");
		return return_error_val();
	}
	#endif
	return (mat->row_start);
}

int get_rows_below_ghost_matrix(const_matrix mat){
	#ifdef _ERROR_CHECKING_ON_
	if(mat==NULL){
		fprintf(stderr,"ERROR mat is NULL in get_rows_below_ghost_matrix\n");
		return return_error_val();
	}
	#endif
	return (mat->rows-mat->row_end);
}
/* start_row, end_row, start_col and end_col can be any
 * value between and including 0 and row-1, col-1 of the
 * smallest of the two matrices */
int copy_matrix_elems(int start_row_elem   ,
                      int end_row_elem     ,
                      int start_col_elem   ,
                      int end_col_elem     ,
                      const_matrix mat_orig,
                      matrix * mat_copy    ){

	int i;
	int j;
	int col_length;
	int row_length;
	float val;

	/* Add 1 because the elements at 0 are included */
	col_length = end_col_elem-start_col_elem+1;
	row_length = end_row_elem-start_row_elem+1;

	#ifdef _ERROR_CHECKING_ON_
	if(mat_copy == NULL){
		fprintf(stderr,"ERROR mat_copy is NULL "
                       "in copy_matrix_elems\n");
		return return_error_val();
	}
	if(mat_orig == NULL){
		fprintf(stderr,"ERROR mat_orig is NULL "
                       "in copy_matrix_elems\n");
		return return_error_val();
	}
	if(start_row_elem>end_row_elem){
		fprintf(stderr,"ERROR start_row_elem is larger than end_row "
                       "in copy_matrix_elems\n");
		return return_error_val();
	}
	if(start_col_elem>end_col_elem){
		fprintf(stderr,"ERROR start_col_elem  larger than end_col in"
                       "copy_matrix_elems\n");
		return return_error_val();
	}
	if(start_row_elem<0){
		fprintf(stderr,"ERROR start_row_elem  is less than 0 in "
                       "copy_matrix_elems\n");
		return return_error_val();
	}
	if(end_row_elem>=mat_orig->rows){
		fprintf(stderr,"ERROR end_row_elem  extends past the row length "
                       "of the original matrix in copy_matrix_elems\n");
		return return_error_val();
	}
	if(row_length>(*mat_copy)->rows){
		fprintf(stderr,"ERROR rows to be copied extend past the "
                       "row length of mat_copy in copy_matrix_elems\n");
		return return_error_val();
	}
	if(start_col_elem<0){
		fprintf(stderr,"ERROR start_col_elem  is less than 0 in "
                       "copy_matrix_elems\n");
		return return_error_val();
	}
	if(end_col_elem>=mat_orig->cols){
		fprintf(stderr,"ERROR end_col_elem  is larger than the column "
                        "length of the originam matrix in "
                        "copy_matrix_elems\n");
		return return_error_val();
	}
	if(col_length>(*mat_copy)->cols){
		fprintf(stderr,"ERROR columns extend past the column"
                       "length of mat_copy in copy_matrix_elems\n");
		return return_error_val();
	}
	#endif

	for(i=0;i<row_length;i++){
		for(j=0;j<col_length;j++){
			val = ELEM(mat_orig,i+start_row_elem,j+start_col_elem);
			ELEM((*mat_copy),i,j) = val;	
		}
	}	

	return 0;
}

float get_elem(const_matrix mat, int row, int col){

	#ifdef _ERROR_CHECKING_ON_
	if(mat==NULL){
		fprintf(stderr,"ERROR mat is NULL in get_elem\n");	
		return return_error_val();
	}
	if(row<0){
		fprintf(stderr,"ERROR row is less than 0 in get_elem\n");
		return return_error_val();
	}
	if(col<0){
		fprintf(stderr,"ERROR col is less than 0 in get_elem\n");
		return return_error_val();
	}
	if(row>=mat->rows){
		fprintf(stderr,"ERROR row greater than rows in matrix get_elem\n");
		return return_error_val();
	}
	if(col>=mat->cols){
		fprintf(stderr,"ERROR col greater than cols in matrix get_elem\n");
		return return_error_val();
	}
	#endif
	return ELEM(mat,row,col);
}

int set_elem(matrix mat, int row, int col, float val){

	#ifdef _ERROR_CHECKING_ON_
	if(row<0){
		fprintf(stderr,"ERROR row is less than 0 in set_elem\n");
		return return_error_val();
	}
	if(col<0){
		fprintf(stderr,"ERROR col is less than 0 in set_elem\n");
		return  return_error_val();
	}
	if(row>=mat->rows){
		fprintf(stderr,"ERROR row greater than rows in matrix set_elem\n");
		return  return_error_val();
	}
	if(col>=mat->cols){
		fprintf(stderr,"ERROR col greater than cols in matrix set_elem\n");
		return  return_error_val();
	}
	#endif
	ELEM(mat,row,col) = val;

	return 0;
}

int set_all_elems(matrix mat, float val){

	#ifdef _ERROR_CHECKING_ON_
	if(mat==NULL){
		fprintf(stderr,"ERROR mat is NULL in set_all_elems\n");
		return return_error_val();
	}
	#endif

	int row;
	int col;

	for(row=0;row<mat->rows;row++){
		for(col=0;col<mat->cols;col++){
			ELEM(mat,row,col) = val;
		}
	}
	return 0;
}

int set_elem_core_matrix(matrix mat, int row, int col, float val){

	#ifdef _ERROR_CHECKING_ON_
	if(mat==NULL){
		fprintf(stderr,"ERROR mat is NULL in set_val_core_matrix\n");
		return return_error_val();
	}
	if(row<0){
		fprintf(stderr,"ERROR row is less than 0 in "
                       "set_elem_core_matrix\n");
		return return_error_val();
	}
	if(col<0){
		fprintf(stderr,"ERROR col is less than 0 in "
                       "set_elem_core_matrix\n");
		return return_error_val();
	}
	if(row>=(mat->row_end-mat->row_start)){
		fprintf(stderr,"ERROR row greater than rows in core matrix "
                       "set_elem_core_matrix\n");
		return return_error_val();
	}
	if(col>=(mat->col_end-mat->col_start)){
		fprintf(stderr,"ERROR col greater than cols in matrix "
                       "set_elem_core_matrix\n");
		return return_error_val();
	}
	#endif
	ELEM(mat,row+mat->row_start,col+mat->col_start) = val;

	return 0;
}

float get_elem_core_matrix(const_matrix mat, int row, int col){

	#ifdef _ERROR_CHECKING_ON_
	if(mat==NULL){
		fprintf(stderr,"ERROR mat is NULL in set_val_core_matrix\n");
		return (float) return_error_val();
	}
	if(row<0){
		fprintf(stderr,"ERROR row is less than 0 in "
                       "get_elem_core_matrix\n");
		return (float) return_error_val();
	}
	if(col<0){
		fprintf(stderr,"ERROR col is less than 0 in "
                       "get_elem_core_matrix\n");
		return (float) return_error_val();
	}
	if(row>=(mat->row_end-mat->row_start)){
		fprintf(stderr,"ERROR row greater than rows in core matrix "
                       "get_elem_core_matrix\n");
		return (float) return_error_val();
	}
	if(col>=(mat->col_end-mat->col_start)){
		fprintf(stderr,"ERROR col greater than cols in matrix "
                       "set_elem_core_matrix\n");
		return (float) return_error_val();
	}
	#endif
	return ELEM(mat,row+mat->row_start,col+mat->col_start);
}

float sum_all_core_matrix_elems(const_matrix mat){

	#ifdef _ERROR_CHECKING_ON_
	if(mat==NULL){
		fprintf(stderr,"ERROR mat is NULL in sum_all_core_matrix\n");
		return return_error_val();
	}
	#endif

	int i;
	int j;
	float sum;
	sum = 0;
	for(i=mat->row_start;i<mat->row_end;i++){
		for(j=mat->col_start;j<mat->col_end;j++){
			sum += ELEM(mat,i,j);
		}
	}

	return sum;
}

int set_elem_top_row_ghost_matrix(matrix mat, int row, int col, float val){

	#ifdef _ERROR_CHECKING_ON_
	if(mat==NULL){
		fprintf(stderr,"ERROR mat is NULL in "
                       "set_val_top_row_ghost_matrix\n");
		return return_error_val();
	}
	if(row<0){
		fprintf(stderr,"ERROR row is less than 0 in "
                       "set_elem_top_row_ghost_matrix\n");
		return return_error_val();
	}
	if(col<0){
		fprintf(stderr,"ERROR col is less than 0 in "
                       "set_elem_top_row_ghost_matrix\n");
		return return_error_val();
	}
	if(mat->row_start==0){
		fprintf(stderr,"ERROR mat does not have a top ghost row in "
                       "set_elem_top_row_ghost_matrix\n");
		return return_error_val();
	}
	if(row>=mat->row_start){
		fprintf(stderr,"ERROR row is larger than the rows in the top"
                       " ghost matrix in get_elem_top_row_ghost_matrix\n");
		return return_error_val();
	}
	if(col>=mat->cols){
		fprintf(stderr,"ERROR col greater than cols in matrix "
                       "set_elem_top_row_ghost_matrix\n");
		return return_error_val();
	}
	#endif
	ELEM(mat,row,col) = val;

	return 0;
}

float get_elem_top_row_ghost_matrix(const_matrix mat, int row, int col){

	#ifdef _ERROR_CHECKING_ON_
	if(mat==NULL){
		fprintf(stderr,"ERROR mat is NULL in "
                       "set_val_top_row_ghost_matrix\n");
		return return_error_val();
	}
	if(row<0){
		fprintf(stderr,"ERROR row is less than 0 in "
                       "get_elem_top_row_ghost_matrix\n");
		return return_error_val();
	}
	if(col<0){
		fprintf(stderr,"ERROR col is less than 0 in "
                       "get_elem_top_row_ghost_matrix\n");
		return return_error_val();
	}
	if(mat->row_start==0){
		fprintf(stderr,"ERROR mat does not have a top ghost row in "
                       "get_elem_top_row_ghost_matrix\n");
		return return_error_val();
	}
	if(row>=mat->row_start){
		fprintf(stderr,"ERROR row is larger than the rows in the top"
                       " ghost matrix in get_elem_top_row_ghost_matrix\n");
		return return_error_val();
	}
	if(col>=mat->cols){
		fprintf(stderr,"ERROR col greater than cols in matrix "
                       "get_elem_top_row_ghost_matrix\n");
		return return_error_val();
	}
	#endif
	return ELEM(mat,row,col);
}

int set_elem_bottom_row_ghost_matrix(matrix mat, int row, int col, float val){

	#ifdef _ERROR_CHECKING_ON_
	if(mat==NULL){
		fprintf(stderr,"ERROR mat is NULL in "
                       "set_val_bottom_row_ghost_matrix\n");
		return return_error_val();
	}
	if(row<0){
		fprintf(stderr,"ERROR row is less than 0 in "
                       "set_elem_bottom_row_ghost_matrix\n");
		return return_error_val();
	}
	if(col<0){
		fprintf(stderr,"ERROR col is less than 0 in "
                       "set_elem_bottom_row_ghost_matrix\n");
		return return_error_val();
	}
	if((mat->rows-mat->row_end)==0){
		fprintf(stderr,"ERROR mat does not have a top ghost row in "
                       "set_elem_bottom_row_ghost_matrix\n");
		return return_error_val();
	}
	if(row>=(mat->rows-mat->row_end)){
		fprintf(stderr,"ERROR row is larger than the size of the "
                       "bottom ghost matrix in "
                       "set_elem_bottom_row_ghost_matrix\n");
		return return_error_val();
	}
	if(col>=mat->cols){
		fprintf(stderr,"ERROR col greater than cols in matrix "
                       "set_elem_bottom_row_ghost_matrix\n");
		return return_error_val();
	}
	#endif
	ELEM(mat,row+mat->row_end,col) = val;

	return 0;
}

float get_elem_bottom_row_ghost_matrix(const_matrix mat, int row, int col){

	#ifdef _ERROR_CHECKING_ON_
	if(mat==NULL){
		fprintf(stderr,"ERROR mat is NULL in "
                       "get_val_bottom_row_ghost_matrix\n");
		return return_error_val();
	}
	if(row<0){
		fprintf(stderr,"ERROR row is less than 0 in "
                       "get_elem_bottom_row_ghost_matrix\n");
		return return_error_val();
	}
	if(col<0){
		fprintf(stderr,"ERROR col is less than 0 in "
                       "get_elem_bottom_row_ghost_matrix\n");
		return return_error_val();
	}
	if((mat->rows-mat->row_end)==0){
		fprintf(stderr,"ERROR mat does not have a top ghost row in "
                       "get_elem_bottom_row_ghost_matrix\n");
		return return_error_val();
	}
	if(row>=(mat->rows-mat->row_end)){
		fprintf(stderr,"ERROR row is larger than the size of the "
                       "bottom ghost matrix in "
                       "get_elem_bottom_row_ghost_matrix\n");
		return return_error_val();
	}
	if(col>=mat->cols){
		fprintf(stderr,"ERROR col greater than cols in matrix "
                       "get_elem_bottom_row_ghost_matrix\n");
		return return_error_val();
	}
	#endif
	return ELEM(mat,row+mat->row_end,col);
}

int send_matrix(matrix mat){

	#ifdef _MPI_H_
	MPI_Datatype Matrix_type;	
	MPI_Aint start_address;
	MPI_Aint address;
	
	MPI_Datatype type[7];
	type[0] = MPI_INT;
	type[1] = MPI_INT;
	type[2] = MPI_INT;
	type[3] = MPI_INT;
	type[4] = MPI_INT;
	type[5] = MPI_INT;
	type[6] = MPI_FLOAT;
	
	MPI_Aint displacements[7];
	displacements[0] = 0;

	MPI_Get_address(&(mat->rows),&start_address);
	MPI_Get_address(&(mat->cols),&address);
	displacements[1] = address-start_address;
	MPI_Get_address(&(mat->col_start),&address);
	displacements[2] = address-start_address;
	MPI_Get_address(&(mat->col_end),&address);
	displacements[3] = address-start_address;
	MPI_Get_address(&(mat->row_start),&address);
	displacements[4] = address-start_address;
	MPI_Get_address(&(mat->row_end),&address);
	displacements[5] = address-start_address;
	MPI_Get_address(&(mat->data),&address);
	displacements[6] = address-start_address;

	int blocklen[7];
	blocklen[0] = 1;
	blocklen[1] = 1;
	blocklen[2] = 1;
	blocklen[3] = 1;
	blocklen[4] = 1;
	blocklen[5] = 1;
	blocklen[6] = get_total_elems(mat);
	
	MPI_Type_create_struct(7, blocklen, displacements,type, &Matrix_type );
	MPI_Type_commit(&Matrix_type);

	MPI_Send(&(mat->rows),1, Matrix_type,dest,0,MPI_COMM_WORLD);

	MPI_Type_free(&Matrix_type);
	#endif
	return 0;
}

int receive_matrix(matrix mat){

	#ifdef _ERROR_CHECKING_ON_
	if(mat==NULL){
		fprintf(stderr,"ERROR receiving matrix to location that "
                       "is NULL. It would be advisable to create"
                       " a matrix in memory first, MPI can then "
                       "write to the location when it recieves "
                       "data structure..");
		return return_error_val();
	}
	#endif
	
	#ifdef _MPI_H_
	MPI_Datatype Matrix_type;	
	MPI_Aint start_address;
	MPI_Aint address;
	
	MPI_Datatype type[7];
	type[0] = MPI_INT;
	type[1] = MPI_INT;
	type[2] = MPI_INT;
	type[3] = MPI_INT;
	type[4] = MPI_INT;
	type[5] = MPI_INT;
	type[6] = MPI_FLOAT;
	
	MPI_Aint displacements[7];
	displacements[0] = 0;

	MPI_Get_address(&(mat->rows),&start_address);
	MPI_Get_address(&(mat->cols),&address);
	displacements[1] = address-start_address;
	MPI_Get_address(&(mat->col_start),&address);
	displacements[2] = address-start_address;
	MPI_Get_address(&(mat->col_end),&address);
	displacements[3] = address-start_address;
	MPI_Get_address(&(mat->row_start),&address);
	displacements[4] = address-start_address;
	MPI_Get_address(&(mat->row_end),&address);
	displacements[5] = address-start_address;
	MPI_Get_address(&(mat->data),&address);
	displacements[6] = address-start_address;

	int blocklen[7];
	blocklen[0] = 1;
	blocklen[1] = 1;
	blocklen[2] = 1;
	blocklen[3] = 1;
	blocklen[4] = 1;
	blocklen[5] = 1;
	blocklen[6] = get_total_elems(mat);
	
	MPI_Type_create_struct(7, blocklen, displacements,type, &Matrix_type );
	MPI_Type_commit(&Matrix_type);

	MPI_Recv(&(mat->rows),1, Matrix_type,source,0,MPI_COMM_WORLD,&status);

	MPI_Type_free(&Matrix_type);
	#endif	
	return 0;
}


/* Will send the col from the core matrix and receive
 * the col to the ghost col, if running checker this 
 * function should be called before the send and receive
 * ghost row function */
bool send_recv_ghost_cols(matrix * new_mat){


	#ifdef _ERROR_CHECKING_ON_
	if(new_mat==NULL){
		fprintf(stderr,"ERROR new_mat is NULL in send_recv_ghost_cols\n");
		return_error_val();
		return false;
	}
	#endif
	int tag = 2;
	float * val_ghost;
	float * val_core;
	int col;
	int row;
	int count;
	int blocklength;
	int stride;
	#ifdef _MPI_H_
	MPI_Datatype Column_type;
	#endif
	/* Will send each column, one at a time. First will 
     * send in the direction indicated in the following
     * diagram:
     *
     *         Rank 0 --> Rank 1 --> Rank 2		 
     *         Rank 3 --> Rank 4 --> Rank 5
     *         Rank 6 --> Rank 7 --> Rank 8 		        
     *       							
     */

	if(verbose){
		printf("Matrix to send and receive columns\n");
		print_matrix_compressed(*new_mat);
	}
	
	dest   = grid_index(my_grid_row,(my_grid_col+1+grid_dim_col)%(grid_dim_col));
	source =  grid_index(my_grid_row,(my_grid_col-1+grid_dim_col)%(grid_dim_col));
          	
	int ghost_cols = 0;
	int ghost_cols_left = 0;
	int ghost_cols_right = 0;

	if((*new_mat)->col_start>0){
		ghost_cols_left = (*new_mat)->col_start;
	}
	if(((*new_mat)->cols-(*new_mat)->col_end)>0){
		ghost_cols_right = (*new_mat)->cols-(*new_mat)->col_end;
	}

	if(ghost_cols_left>ghost_cols_right){
		ghost_cols = ghost_cols_left;
	}else{
		ghost_cols = ghost_cols_right;
	}

	/* We must create a derived data type to deal with the columns 
	 * count       - for a column is the number of elements in a 
	 *               column this is equivalent to the number of rows.
	 * blocklength - Are only taking a single element at each row in
	 *               in the column thus block length will be 1
	 * stride      - The stride is equivalent to the number of rows 
	 *               in the original matrix not the core matrix */
	count       = get_rows_core_matrix(*new_mat);
	blocklength = ghost_cols;
	stride      = get_cols(*new_mat);

	if(verbose){
		printf("Vector type:\n");
		printf("count %d\n",count);
		printf("blocklength %d\n",blocklength);
		printf("stride %d\n",stride);
		printf("Rows %d Cols %d in matrix\n",
				get_rows(*new_mat),get_cols(*new_mat));
	}
	#ifdef _MPI_H_
	MPI_Type_vector(count,blocklength, 
					stride,MPI_FLOAT, 
					&Column_type);

	MPI_Type_commit(&Column_type);
	#endif

	if(my_grid_col==0){
		/* Grab the float located at the last column of
		 * the core new_matrix. We will grab the columns at the 
		 * right of the core new_matrix because they are
		 * what we will be sending */
		col = (*new_mat)->col_end-ghost_cols;
		row = (*new_mat)->row_start;
		if(verbose){
			printf("Rank %d sending columns starting at ELEM(new_mat,%d,%d)=%f"
					" sending to %d\nNote if row is not 0 it just means "
					"that there are ghost rows. When sending the columns"
					" we are excluding the ghost rows thus there may be "
					"an offset which would explain why rows may not be 0"
					"\n",my_rank,row,col,ELEM(*new_mat,row,col),dest); 	
		}

		val_core = &(ELEM(*new_mat,row,col));

		if(send_recv_sync){
			#ifdef _MPI_H_
			MPI_Send(val_core, 1, Column_type, dest, tag,
					MPI_COMM_WORLD);	
			#endif
		}else{
			#ifdef _MPI_H_
			MPI_Isend(val_core, 1, Column_type, dest, tag,
					MPI_COMM_WORLD,&send_request3);	
			#endif
		}
	}else if(my_grid_col==(grid_dim_col-1)){
		/* Grab the float located at the first row 
		 * (in the core matrix) in the left ghost col. 
		 * This is where we will be receiving from */
		col = 0;
		row = (*new_mat)->row_start;
		if(verbose){
			printf("Rank %d receiving value at ELEM(new_mat,%d,%d)=%f from %d\n"
					"Note rows may not be NULL especially if there are ghost"
					" rows above and below the core matrix. When we send the" 
					" columns we are excluding the elements that exist in the"
					" ghost rows.",my_rank,row,col,ELEM(*new_mat,row,col),source);
		}
		val_ghost=&(ELEM((*new_mat),row,col));

		if(send_recv_sync){
			#ifdef _MPI_H_
			MPI_Recv(val_ghost, 1, Column_type, source, tag, 
					 MPI_COMM_WORLD, &status);	
			#endif
		}else{
			#ifdef _MPI_H_
			MPI_Irecv(val_ghost, 1, Column_type, source, tag, 
					  MPI_COMM_WORLD, &recv_request3);	
			#endif
		}

	}else{
		/* Grab the float located at the last column of
		 * the core new_matrix. We will grab the columns at the 
		 * right of the core new_matrix because they are
		 * what we will be sending */
		col = (*new_mat)->col_end-ghost_cols;
		row = (*new_mat)->row_start;
		val_core = &(ELEM((*new_mat),row,col));
		/* Grab the float located at the first row 
		 * (in the core matrix) in the left ghost col. 
		 * This is where we will be receiving from */
		col = 0;
		row = (*new_mat)->row_start;
		val_ghost=&(ELEM((*new_mat),row,col));

		if(verbose){
			printf("Rank %d Sending to %d and Receiving from %d\n",
					my_rank, dest, source);
		}
		if(send_recv_sync){
			#ifdef _MPI_H_
			MPI_Send(val_core, 1, Column_type, dest, tag,
					 MPI_COMM_WORLD);	
			MPI_Recv(val_ghost, 1, Column_type, source, tag, 
					 MPI_COMM_WORLD, &status);	
			#endif
		}else{
			#ifdef _MPI_H_
			MPI_Isend(val_core, 1, Column_type, dest, tag,
					  MPI_COMM_WORLD, &send_request3);	
			MPI_Irecv(val_ghost, 1, Column_type, source, tag, 
					  MPI_COMM_WORLD, &recv_request3);	
			#endif
		}

	}

	/* Will send each col, one at a time. Will now
	 * send in the direction indicated in the following
	 * diagram:
	 *         Rank 0 <-- Rank 1 <-- Rank 2		 
	 *         Rank 3 <-- Rank 4 <-- Rank 5
	 *         Rank 6 <-- Rank 7 <-- Rank 8 		        
	 *       							
	 */

	/* Switch the source and destination processors
	 * will now be receiving from */
	dest   = grid_index(my_grid_row,(my_grid_col-1+grid_dim_col)%(grid_dim_col));
	source =  grid_index(my_grid_row,(my_grid_col+1+grid_dim_col)%(grid_dim_col));

	/* Increment until we have sent all of the ghost_rows
	 * located on the top of the core new_matrix */
	tag = 3;

	/* Sending left column and receiving to right column */
	if(my_grid_col==0){
		/* Grab the float located at the last column in
		 * the ghost columns and within the core matrix row.
		 * This is where we will be receiving from.*/
		row = (*new_mat)->row_start;
		col = (*new_mat)->col_end;
		if(verbose){
			printf("Rank %d receiving value at ELEM(new_mat,%d,%d)=%f from %d\n",my_rank,row,col,get_elem(*new_mat,row,col),source);
		}

		val_ghost= &(ELEM((*new_mat),row,col));
		if(send_recv_sync){
			#ifdef _MPI_H_
			MPI_Recv(val_ghost, 1, Column_type, source, tag, 
					 MPI_COMM_WORLD, &status);
			#endif
		}else{
			#ifdef _MPI_H_
			MPI_Irecv(val_ghost, 1, Column_type, source, tag, 
					  MPI_COMM_WORLD, &recv_request4);
			#endif
		}
	}else if(my_grid_col==(grid_dim_col-1)){
		/* Grab the float located at the first row of
		 * the core new_matrix and the first column of the core_ 
		 * matrix. We will grab the cols at the left
		 * of the core new_matrix because these are what we will be
		 * sending to the right of the destinations processors ghost
		 * col */
		row = (*new_mat)->row_start;
		col = (*new_mat)->col_start;
		if(verbose){
			printf("Rank %d sending value at ELEM(new_mat,%d,%d)=%f to %d\n",my_rank,row,col,get_elem(*new_mat,row,col),dest);
		}

		val_core=&(ELEM((*new_mat),row,col));
		if(send_recv_sync){
			#ifdef _MPI_H_
			MPI_Send(val_core, 1, Column_type, dest, tag,
					 MPI_COMM_WORLD);	
			#endif
		}else{
			#ifdef _MPI_H_
			MPI_Isend(val_core, 1, Column_type, dest, tag,
					  MPI_COMM_WORLD, &send_request4);	
			#endif
		}
	}else{
		/* Grab the float located at the first row of
		 * the core new_matrix and the first column of the core_ 
		 * matrix. We will grab the cols at the left
		 * of the core new_matrix because these are what we will be
		 * sending to the right of the destinations processors ghost
		 * col */
		row = (*new_mat)->row_start;
		col = (*new_mat)->col_start;
		val_core=&(ELEM((*new_mat),row,col));
		/* Grab the float located at the last column in
		 * the ghost columns and within the core matrix row.
		 * This is where we will be receiving from.*/
		row = (*new_mat)->row_start;
		col = (*new_mat)->col_end;
		val_ghost=&(ELEM((*new_mat),row,col));
		/* Notice that rank 0 does not participate because
		 * it does not have a top ghost row and our World is
		 * not periodic. */
		if(verbose){
			printf("Rank %d Sending to %d and Receiving from %d\n",
					my_rank, dest, source);
		}
		if(send_recv_sync){
			#ifdef _MPI_H_
			MPI_Send(val_core, 1, Column_type, dest, tag,
					 MPI_COMM_WORLD);	
			MPI_Recv(val_ghost, 1, Column_type, source, tag, 
					 MPI_COMM_WORLD, &status);	
			#endif
		}else{
			#ifdef _MPI_H_
			MPI_Isend(val_core, 1, Column_type, dest, tag,
					  MPI_COMM_WORLD, &send_request4);	
			MPI_Irecv(val_ghost, 1, Column_type, source, tag, 
					  MPI_COMM_WORLD, &recv_request4);	
			#endif
		}
	}

	/* Free the column data type */
	#ifdef _MPI_H_
	MPI_Type_free(&Column_type);
	#endif
	if(verbose){
		printf("Matrix after sending and receiving columns\n");
		print_matrix_compressed(*new_mat);
	}
	return true;
}




/* Will send the row from the core matrix and receive
 * the row to the ghost row */
bool send_recv_ghost_rows(matrix * new_mat){

	#ifdef _ERROR_CHECKING_ON_
	if(new_mat==NULL){
		fprintf(stderr,"ERROR new_mat is NULL in send_recv_ghost_rows\n");
		return_error_val();
		return false;
	}
	#endif
	int tag = 0;
	float * val_ghost;
	float * val_core;
	int row;
	/* Will send each row, one at a time. First will 
	 * send in the direction indicated in the following
	 * diagram:
	 *
	 *      |  Rank 0				top panel (panel 1) of matrix 
	 *      V | Rank 1				          (panel 2) of matrix
	 *        V | Rank 2							:
	 *          V | ...								:
	 *            V Rank (num_proc-1)		  (panel num_proc) of matrix
	 */

	if(verbose){
		printf("Matrix to send and receive rows\n");
		print_matrix_compressed(*new_mat);
	}
	dest   = grid_index((my_grid_row+1+grid_dim_row)%(grid_dim_row),my_grid_col);
	source = grid_index((my_grid_row-1+grid_dim_row)%(grid_dim_row),my_grid_col);

	int ghost_rows = 0;
	int ghost_rows_top = 0;
	int ghost_rows_bottom = 0;

	if((*new_mat)->row_start>0){
		ghost_rows_top = (*new_mat)->row_start;
	}
	if(((*new_mat)->rows-(*new_mat)->row_end)>0){
		ghost_rows_bottom = (*new_mat)->rows-(*new_mat)->row_end;
	}

	if(ghost_rows_top>ghost_rows_bottom){
		ghost_rows = ghost_rows_top;
	}else{
		ghost_rows = ghost_rows_bottom;
	}

	if(my_grid_row==0){
		/* Grab the float located at the first column of
		 * the core new_matrix. We will grab the rows at the 
		 * bottom of the core new_matrix because there are
		 * what we will be sending */
		row = (*new_mat)->row_end-ghost_rows;
		if(verbose){
			printf("Rank %d sending value at ELEM(new_mat,%d,0)=%f to %d\n",my_rank,row, ELEM(*new_mat,row,0),dest);
			printf("Sending total of %d columns\n",(*new_mat)->cols);
		}	

		val_core = &(ELEM(*new_mat,row,0));
		if(send_recv_sync){
			#ifdef _MPI_H_
			MPI_Send(val_core, (*new_mat)->cols*ghost_rows, MPI_FLOAT, dest, tag,
					MPI_COMM_WORLD);	
			#endif
		}else{
			#ifdef _MPI_H_
			MPI_Isend(val_core, (*new_mat)->cols*ghost_rows, MPI_FLOAT, dest, tag,
					MPI_COMM_WORLD,&send_request1);	
			#endif
		}
	}else if(my_grid_row==(grid_dim_row-1)){
		/* Grab the float located at the first column in
		 * the bottom ghost row. This is where we will be 
		 * receiving from */
		row = 0;
		if(verbose){
			printf("Rank %d receiving value at ELEM(new_mat,%d,0)=%f from %d\n",my_rank,row,ELEM(*new_mat,row,0),source);
			printf("Sending total of %d columns\n",(*new_mat)->cols);
		}

		val_ghost=&(ELEM((*new_mat),row,0));
		if(send_recv_sync){
			#ifdef _MPI_H_
			MPI_Recv(val_ghost, (*new_mat)->cols*ghost_rows, MPI_FLOAT, source, tag, 
					MPI_COMM_WORLD, &status);	
			#endif
		}else{
			#ifdef _MPI_H_
			MPI_Irecv(val_ghost, (*new_mat)->cols*ghost_rows, MPI_FLOAT, source, tag, 
					MPI_COMM_WORLD, &recv_request1);	
			#endif
		}

	}else{
		/* Grab the float located at the first column of
		 * the core new_matrix. We will grab the rows at the 
		 * bottom of the core new_matrix because there are
		 * what we will be sending */
		row = (*new_mat)->row_end-ghost_rows;
		val_core = &(ELEM((*new_mat),row,0));
		/* Grab the float located at the first column in
		 * the bottom ghost row. This is where we will be 
		 * receiving from */
		row = 0;
		val_ghost=&(ELEM((*new_mat),row,0));
		/* If rank 0 sending and receiving to proc 1 
		 * sending bottom row of core new_matrix of proc 0 
		 * to top ghost row of proc 1. Note that the 
		 * last rank does not participate */
		if(verbose){
			printf("Rank %d Sending to %d and Receiving from %d\n",
					my_rank, dest, source);
			printf("Sending total of %d columns\n",(*new_mat)->cols);
		}
		if(send_recv_sync){
			#ifdef _MPI_H_
			MPI_Send(val_core, (*new_mat)->cols*ghost_rows, MPI_FLOAT, dest, tag,
					MPI_COMM_WORLD);	
			MPI_Recv(val_ghost, (*new_mat)->cols*ghost_rows, MPI_FLOAT, source, tag, 
					MPI_COMM_WORLD, &status);	
			#endif
		}else{
			#ifdef _MPI_H_
			MPI_Isend(val_core, (*new_mat)->cols*ghost_rows, MPI_FLOAT, dest, tag,
					MPI_COMM_WORLD, &send_request1);	
			MPI_Irecv(val_ghost, (*new_mat)->cols*ghost_rows, MPI_FLOAT, source, tag, 
					MPI_COMM_WORLD, &recv_request1);	
			#endif
		}

	}

	/* Will send each row, one at a time. Will now
	 * send in the direction indicated in the following
	 * diagram:
	 *
	 *      A  Rank 0				top panel (panel 1) of new_matrix 
	 *      | A Rank 1				          (panel 2) of new_matrix
	 *        | A Rank 2							:
	 *          | A ...								:
	 *            | Rank (num_proc-1)		  (panel num_proc) of new_matrix
	 */

	/* Switch the source and destination processors
	 * will now be receiving from */
	dest   = grid_index((my_grid_row-1+grid_dim_row)%(grid_dim_row),my_grid_col);
	source = grid_index((my_grid_row+1+grid_dim_row)%(grid_dim_row),my_grid_col);

	/* Increment until we have sent all of the ghost_rows
	 * located on the top of the core new_matrix */
	tag = 1;

	/* Sending top row and receiving to bottom row */
	if(my_grid_row==0){
		/* Grab the float located at the first column in 
		 * the top ghost row. This is where we will be receiving
		 * from.*/
		/* Note we minus 1 from the ghost row here because we are
		 * starting on the ghost row with (*new_mat)->row_end*/
		row = (*new_mat)->row_end;
		if(verbose){
			printf("Rank %d receiving value at ELEM(new_mat,%d,0)=%f from %d\n",my_rank,row,get_elem(*new_mat,row,0),source);
			printf("Sending total of %d columns\n",(*new_mat)->cols);
		}
		val_ghost= &(ELEM((*new_mat),row,0));
		if(send_recv_sync){
			#ifdef _MPI_H_
			MPI_Recv(val_ghost, (*new_mat)->cols*ghost_rows, MPI_FLOAT, source, tag, 
					MPI_COMM_WORLD, &status);
			#endif
		}else{
			#ifdef _MPI_H_
			MPI_Irecv(val_ghost, (*new_mat)->cols*ghost_rows, MPI_FLOAT, source, tag, 
					MPI_COMM_WORLD, &recv_request2);
			#endif
		}
	}else if(my_grid_row==(grid_dim_row-1)){
		/* Grab the float located at the first column of
		 * the core new_matrix. We will grab the rows at the bottom
		 * of the core new_matrix because these are what we will be
		 * sending to the top of the destinations processors ghost
		 * row */
		row = (*new_mat)->row_start;
		if(verbose){
			printf("Rank %d sending value at ELEM(new_mat,%d,0)=%f to %d\n",my_rank,row,get_elem(*new_mat,row,0),dest);
			printf("Sending total of %d columns\n",(*new_mat)->cols);
		}
		val_core=&(ELEM((*new_mat),row,0));
		if(send_recv_sync){
			#ifdef _MPI_H_
			MPI_Send(val_core, (*new_mat)->cols*ghost_rows, MPI_FLOAT, dest, tag,
					MPI_COMM_WORLD);	
			#endif
		}else{
			#ifdef _MPI_H_
			MPI_Isend(val_core, (*new_mat)->cols*ghost_rows, MPI_FLOAT, dest, tag,
					MPI_COMM_WORLD, &send_request2);	
			#endif
		}
	}else{
		/* Grab the float located at the first column of
		 * the core new_matrix. We will grab the rows at the bottom
		 * of the core new_matrix because these are what we will be
		 * sending to the top of the destinations processors ghost
		 * row */
		row = (*new_mat)->row_start;
		val_core=&(ELEM((*new_mat),row,0));
		/* Grab the float located at the first column in 
		 * the top ghost row. This is where we will be receiving
		 * from.*/
		row = (*new_mat)->row_end;
		val_ghost=&(ELEM((*new_mat),row,0));
		/* Notice that rank 0 does not participate because
		 * it does not have a top ghost row and our World is
		 * not periodic. */
		if(verbose){
			printf("Rank %d Sending to %d and Receiving from %d\n",
					my_rank, dest, source);
			printf("Sending total of %d columns\n",(*new_mat)->cols);
		}
		if(send_recv_sync){
			#ifdef _MPI_H_
			MPI_Send(val_core, (*new_mat)->cols*ghost_rows, MPI_FLOAT, dest, tag,
					MPI_COMM_WORLD);	
			MPI_Recv(val_ghost, (*new_mat)->cols*ghost_rows, MPI_FLOAT, source, tag, 
					MPI_COMM_WORLD, &status);	
			#endif
		}else{
			#ifdef _MPI_H_
			MPI_Isend(val_core, (*new_mat)->cols*ghost_rows, MPI_FLOAT, dest, tag,
					MPI_COMM_WORLD, &send_request2);	
			MPI_Irecv(val_ghost, (*new_mat)->cols*ghost_rows, MPI_FLOAT, source, tag, 
					MPI_COMM_WORLD, &recv_request2);	
			#endif
		}
	}

	if(verbose){
		printf("Matrix after sending and receiving rows\n");
		print_matrix_compressed(*new_mat);
	}
	return true;
}



float live_cell(float sum, float cell){

	#ifdef _ERROR_CHECKING_ON_
	if(sum<0.0){
		fprintf(stderr,"ERROR sum is less than 0.0 should be "
                       "positive or 0.0 in live_cell\n");
		return (float) return_error_val();
	}
	if(sum>8.0){
		fprintf(stderr,"ERROR sum is greater than 8.0, should "
                       "at max be 8.0, because there are only "
                       "8.0 neighbors in live_cell\n");
		return (float) return_error_val();
	}
	if(cell<0.0){
		fprintf(stderr,"ERROR cell is negative can only be 0 or"
                       " 1 in live_cell\n");
		return (float) return_error_val();
	}
	if(cell>1.0){
		fprintf(stderr,"ERROR cell is greater than 1.0, can "
                       "only be 0 or 1 in live_cell\n");
		return (float) return_error_val();
	}
	#endif

	switch((int)sum){
		case 2:
			if(cell==1.0){
				return 1.0;
			}else{
				return 0.0;
			}
		case 3:
			return 1.0;
	}

	return 0.0;	
}

static inline void conways_condition(const_matrix orig_mat, 
                                     matrix * new_mat, 
                                     int i, int j,
                                     int ghost){

	#ifdef _ERROR_CHECKING_ON_
	if(orig_mat==NULL){
		fprintf(stderr,"ERROR orig_mat is NULL in conways_condition\n");
		return_error_val();
	}
	if(new_mat==NULL){
		fprintf(stderr,"ERROR new_mat is NULL in conways_condition\n");
		return_error_val();
	}
	if((orig_mat)->rows!=(*new_mat)->rows){
		fprintf(stderr,"ERROR orig_mat and new_mat do not have the same"
                       " number of rows in conways_condition\n");
		fprintf(stderr,"Rank %d\n",my_rank);
		fprintf(stderr,"orig_mat %d rows\n",(orig_mat)->rows);
		fprintf(stderr,"new_mat  %d rows\n",(*new_mat)->rows);
		return_error_val();
	}
	if((orig_mat)->cols!=(*new_mat)->cols){
		fprintf(stderr,"ERROR orig_mat and new_mat do not have the same"
                       " number of columns in conways_condition\n");
		return_error_val();
	}
	if((orig_mat)->row_start!=(*new_mat)->row_start){
		fprintf(stderr,"ERROR orig_mat and new_mat do not have the same"
                       " number of top ghost rows in conways_condition\n");
		return_error_val();
	}
	if((orig_mat)->row_end!=(*new_mat)->row_end){
		fprintf(stderr,"ERROR orig_mat and new_mat do not have the same "
                       "number bottom ghost rows in conways_condition\n");
		return_error_val();
	}
	if(block_on){
		if((*new_mat)->col_start!=0){
			fprintf(stderr,"ERROR new_mat has left ghost columns which "
                           "are not needed in the panel implementation "
                           "(block) in conways_condition\n");
			return_error_val();
		}
		if((*new_mat)->col_end!=(*new_mat)->cols){
			fprintf(stderr,"ERROR new_mat has right ghost columns which "
                           "are not needed in the panel implementation "
                           "(block) in conways_condition\n");
			return_error_val();
		}	
		if((orig_mat)->col_start!=0){
			fprintf(stderr,"ERROR orig_mat has left ghost columns which "
                           "are not needed in the panel implementation "
                           "(block) in conways_condition\n");
			return_error_val();
		}
		if((orig_mat)->col_end!=(orig_mat)->cols){
			fprintf(stderr,"ERROR orig_mat has right ghost columns which "
                           "are not needed in the panel implementation "
                           "(block) in conways_condition\n");
			return_error_val();
		}
	}
	if(i<0){
		fprintf(stderr,"ERROR i is less than 0 in conways_condition\n");
		return_error_val();
	}
	if(j<0){
		fprintf(stderr,"ERROR j is less than 0 in conways_condition\n");
		return_error_val();
	}
	if(i>=(orig_mat)->rows){
		fprintf(stderr,"ERROR i is larger than or equal to the number of "
                       "rows in *orig_mat in conways_condition\n");
		return_error_val();
	}
	if(j>=(orig_mat)->cols){
		fprintf(stderr,"ERROR j is larger than or equal to the number of "
                       "cols in *orig_mat in conways_condition\n");
		return_error_val();
	}
	#endif
	
	float value;

	float NW;
	float N;
	float NE;
	float E;
	float W;
	float SW;
	float S;
	float SE;

	float sum = 0.0;
	value = get_elem(orig_mat,i,j);
	if(within_matrix(orig_mat,i-1,j-1)){
		NW = get_elem(orig_mat,i-1,j-1);
		sum += (int)NW;
	}
	if(within_matrix(orig_mat,i,j-1)){
		N = get_elem(orig_mat,i,j-1);
		sum += (int)N;
	}
	if(within_matrix(orig_mat,i+1,j-1)){
		NE = get_elem(orig_mat,i+1,j-1);
		sum += (int)NE;
	}
	if(within_matrix(orig_mat,i-1,j)){
		W = get_elem(orig_mat,i-1,j);
		sum += (int)W;
	}
	if(within_matrix(orig_mat,i+1,j)){
		E = get_elem(orig_mat,i+1,j);
		sum += (int)E;
	}
	if(within_matrix(orig_mat,i-1,j+1)){
		SW = get_elem(orig_mat,i-1,j+1);
		sum += (int)SW;
	}
	if(within_matrix(orig_mat,i,j+1)){
		S = get_elem(orig_mat,i,j+1);
		sum += (int)S;
	}
	if(within_matrix(orig_mat,i+1,j+1)){
		SE = get_elem(orig_mat,i+1,j+1);
		sum += (int)SE;
	}

	if(within_core_matrix_with_ghost((*new_mat),i,j,ghost)){
		set_elem(*new_mat,i,j,live_cell(sum, value));
	}

}

bool within_core_matrix(const_matrix mat, int i, int j){

	#ifdef _ERROR_CHECKING_ON_
	if(mat==NULL){
		fprintf(stderr,"ERROR mat is NULL in within_core_matrix\n");
		return_error_val();
		return false;
	}
	if(i<0){
		fprintf(stderr,"ERROR i is less than 0 in within_core_matrix\n");
		return_error_val();
		return false;
	}
	if(j<0){
		fprintf(stderr,"ERROR j is less than 0 in within_core_matrix\n");
		return_error_val();
		return false;
	}
	if(i>=mat->rows){
		fprintf(stderr,"ERROR i is larger than or equal to the number of"
                       " rows in mat in within_core_matrix\n");
		return_error_val();
		return false;
	}
	if(j>=mat->cols){
		fprintf(stderr,"ERROR j is larger than or equal to the number of"
                       " cols in mat in within_core_matrix\n");
		return_error_val();
		return false;
	}
	#endif

	if(i<mat->row_start){
		return false;
	}else if(i>=mat->row_end){
		return false;
	}else if(j<mat->col_start){
		return false;
	}else if(j>=mat->col_end){
		return false;
	}else{
		return true;
	}

}

bool within_core_matrix_with_ghost(const_matrix mat, int i, int j, int ghost){

	#ifdef _ERROR_CHECKING_ON_
	if(mat==NULL){
		fprintf(stderr,"ERROR mat is NULL in within_core_matrix_with_ghost\n");
		return_error_val();
		return false;
	}
	if(i<0){
		fprintf(stderr,"ERROR i is less than 0 in within_core_matrix_with_ghost\n");
		return_error_val();
		return false;
	}
	if(j<0){
		fprintf(stderr,"ERROR j is less than 0 in within_core_matrix_with_ghost\n");
		return_error_val();
		return false;
	}
	if(i>=mat->rows){
		fprintf(stderr,"ERROR i is larger than or equal to the number of"
                       " rows in mat in within_core_matrix_with_ghost\n");
		return_error_val();
		return false;
	}
	if(j>=mat->cols){
		fprintf(stderr,"ERROR j is larger than or equal to the number of"
                       " cols in mat in within_core_matrix_with_ghost\n");
		return_error_val();
		return false;
	}
	if(ghost<0){
		fprintf(stderr,"ERROR ghost is less than 0 in "
                       "within_core_matrix_with_ghost\n");
		return_error_val();
		return false;
	}
	#endif

	if(i<(mat->row_start-ghost)){
		return false;
	}else if(i>=(mat->row_end+ghost)){
		return false;
	}else if(j<(mat->col_start-ghost)){
		return false;
	}else if(j>=(mat->col_end+ghost)){
		return false;
	}else{
		return true;
	}

}


bool within_matrix(const_matrix mat, int i, int j){

	#ifdef _ERROR_CHECKING_ON_
	if(mat==NULL){
		fprintf(stderr,"ERROR mat is NULL in within_core_matrix\n");
		return_error_val();
		return false;
	}
	#endif
	
	if(i<0){
		return false;
	}else if(j<0){
		return false;
	}else if(i>=mat->rows){
		return false;
	}else if(j>=mat->cols){
		return false;
	}else{
		return true;
	}
}
int update_border_conway(matrix temp1 , matrix * temp2, int k){

	#ifdef _ERROR_CHECKING_ON_
	if(temp1==NULL){
		fprintf(stderr,"ERROR temp1 is NULL in "
                       "update_border_conway\n");
		return return_error_val();
	}
	if(temp2==NULL){
		fprintf(stderr,"ERROR temp2 is NULL in "
                       "update_border_conway\n");
		return return_error_val();
	}
	if(get_rows(temp1)!=get_rows(*temp2)){
		fprintf(stderr,"ERROR rows of temp1 %d and temp2 %d are "
                       "not equal in update_border_conway\n",
                        get_rows(temp1),get_rows(*temp2));
		return return_error_val();
	}
	if(get_cols(temp1)!=get_cols(*temp2)){
		fprintf(stderr,"ERROR cols of temp1 and temp2 are "
                       "not equal in update_border_conway\n");
		return return_error_val();
	}
	if(k<0){
		fprintf(stderr,"ERROR k is less than 0 in "
                       "update_border_conway\n");
		return return_error_val();
	}
	#endif

	int row_end;
	int row_start;
	int col_end;
	int col_start;

	int ghost_left  = get_cols_left_ghost_matrix(temp1);
	int ghost_right = get_cols_right_ghost_matrix(temp1);
	int ghost_above = get_rows_above_ghost_matrix(temp1);
	int ghost_below = get_rows_below_ghost_matrix(temp1);
	
	/* Updating Left Border West */
	if(ghost_left>0){
		if(ghost_below){
			row_end = get_row_end(temp1)-(ghost_below)*2+1+k;
		}else{
			row_end = get_rows(temp1);
		}		
		if(ghost_above){
			row_start = get_row_start(temp1)+(ghost_above*2-1)-k;
		}else{
			row_start = 0;
		}
		col_start = k+1;
		col_end   = get_col_start(temp1)+ghost_left*2-1-k;	
		if(verbose){
			printf("Updating Left Border West\n");
			printf("col_start %d col_end %d\n",col_start,col_end);
			printf("row_start %d row_end %d\n",row_start,row_end);
			printf("ghost_left %d k %d get_col_start %d\n",ghost_left,k,get_col_start(temp1));
			fflush(NULL);
		}
		set_elems_conway_rectangle(col_start,col_end,
								   row_start,row_end,
								   temp1, temp2,ghost_left-k);
		if(verbose){
			print_matrix_compressed(*temp2);
			fflush(NULL);
		}
	}

	/* Updating Right Border East */
	if(ghost_right>0){
		if(ghost_below){
			row_end = get_row_end(temp1)-(ghost_below)*2+1+k;
		}else{
			row_end = get_rows(temp1);
		}
		if(ghost_above){
			row_start = get_row_start(temp1)-(ghost_above-1)+k;
		}else{
			row_start = 0;
		}
		col_start = get_col_end(temp1)-ghost_right*2+k+1;
		
		col_end   = get_col_end(temp1)+ghost_right-1-k;
		if(verbose){
			printf("Updating Right Border East\n");
			printf("col_start %d col_end %d\n",col_start,col_end);
			printf("row_start %d row_end %d\n",row_start,row_end);
			printf("ghost_right %d k %d get_col_start %d\n",ghost_right,k,get_col_start(temp1));
			fflush(NULL);
		}
		set_elems_conway_rectangle(col_start,col_end,
								   row_start,row_end,
								   temp1, temp2,ghost_right-k);
		if(verbose){
			print_matrix_compressed(*temp2);
			fflush(NULL);
		}
	}
	/* Top Row North */
	if(ghost_above){
		/* Outer Margin */
		if(ghost_left==0){
			col_start=0;
		}else{
			col_start = 1+k;
		}
		if(ghost_right==0){
			col_end = get_cols(temp1);
		}else{
			col_end = get_cols(temp1)-k-1;
		}
		row_start = 1+k;
		row_end   = get_row_start(temp1)+(ghost_above*2)-k-1;
		set_elems_conway_rectangle(col_start,col_end,
								   row_start,row_end,
								   temp1, temp2, 
                                   ghost_above-k);
		if(verbose){
			printf("Updating North Row Outer Margin\n");
			print_matrix_compressed(*temp2);
			fflush(NULL);
		}
	}
	/* Bottom Row South */ 
	if(ghost_below){
		/* Outer Margin */
		if(ghost_left==0){
			col_start = 0;
		}else{
			col_start = 1+k;
		}
		if(ghost_right==0){
			col_end = get_cols(temp1);
		}else{
			col_end = get_cols(temp1)-k-1;
		}

		row_start = get_row_end(temp1)-(ghost_below*2-k-1);
		row_end   = get_row_end(temp1)+(ghost_below-k-2)+1;
		if(verbose){
			printf("Bottom Row South Outer Margin\n");
			printf("col_start %d col_end %d\n",col_start,col_end);
			printf("row_start %d row_end %d\n",row_start,row_end);
			printf("get_end %d ghost_below %d k %d\n",get_row_end(temp1),ghost_below,k);
			fflush(NULL);
		}
		set_elems_conway_rectangle(col_start,col_end,
				                   row_start,row_end,
				                   temp1, temp2,
                                   ghost_below-k);	
		if(verbose){
			print_matrix_compressed(*temp2);
			fflush(NULL);
		}
	}

	return 0;
}

int generate_data_file(char * filename, matrix mat){

	int i;
	int j;
	FILE * fid;
	fid = fopen(filename,"w");

	for(j=0;j<get_cols(mat);j++){
		for(i=0;i<get_rows(mat);i++){
			fprintf(fid,"%6d ",(int)get_elem(mat,i,j));
		}
		fprintf(fid,"\n");
	}
	fclose(fid);
	return 0;
}

int update_conways_right(matrix orig, matrix temp, matrix new_mat, int k){

	#ifdef _ERROR_CHECKING_ON_
	if(orig==NULL){
		fprintf(stderr,"ERROR orig is NULL in "
                       "update_conways_right_even\n");
		return return_error_val();
	}
	if(temp==NULL){
		fprintf(stderr,"ERROR temp is NULL in "
                       "update_conways_right_even\n");
		return return_error_val();
	}
	if(new_mat==NULL){
		fprintf(stderr,"ERROR new is NULL in "
                       "update_conways_right_even\n");
		return return_error_val();
	}
	if(k<0){
		fprintf(stderr,"ERROR k is less than 0 in "
                       "update_conways_right_even\n");
		return return_error_val();
	}
	if(get_rows(orig)!=get_rows(temp)){
		fprintf(stderr,"ERROR rows in orig is not equal to "
                       "rows in temp in "
                       "update_conways_right_even\n");
		return return_error_val();
	}
	if(get_rows(new_mat)!=get_rows(temp)){
		fprintf(stderr,"ERROR rows in new_mat is not equal to"
                       " rows in temp in "
                       "update_conways_right_even\n");
		return return_error_val();
	}
	if(get_cols(orig)!=get_cols(temp)){
		fprintf(stderr,"ERROR cols in orig is not equal to "
                       "cols in temp in "
                       "update_conways_right_even\n");
		return return_error_val();
	}
	if(get_cols(new_mat)!=get_cols(temp)){
		fprintf(stderr,"ERROR cols in new_mat is not equal "
                       "to columns in temp in "
                       "update_conways_right_even\n");
		return return_error_val();
	}
	#endif

	int row_end;
	int row_start;
	int col_end;
	int col_start;

	int ghost_left  = get_cols_left_ghost_matrix(orig);
	int ghost_right = get_cols_right_ghost_matrix(orig);
	int ghost_above = get_rows_above_ghost_matrix(orig);
	int ghost_below = get_rows_below_ghost_matrix(orig);

	if(ghost_above==0){
		row_start = 0;
	}else{
		row_start = get_row_start(orig)+1+k;
	}
	if(ghost_below==0){
		row_end = get_rows(orig);
	}else{
		row_end = get_row_end(orig)-1-k;
	}
	if(ghost_left==0){
		col_start=0;
	}else{
		col_start = get_col_start(orig)+1+k;
	}
	if(ghost_right==0){
		col_end = get_cols(orig);
	}else{
		col_end = get_col_end(orig)-1-k;
	}

	set_elems_conway_rectangle_rotation(col_end/2, col_end,
										row_start, row_end,
										orig , &temp, &new_mat, k);
	return 0;
}

int update_conways_left(matrix orig, matrix temp, matrix new_mat, int k){

	#ifdef _ERROR_CHECKING_ON_
	if(orig==NULL){
		fprintf(stderr,"ERROR orig is NULL in "
                       "update_conways_left_even\n");
		return return_error_val();
	}
	if(temp==NULL){
		fprintf(stderr,"ERROR temp is NULL in "
                       "update_conways_left_even\n");
		return return_error_val();
	}
	if(new_mat==NULL){
		fprintf(stderr,"ERROR new is NULL in "
                       "update_conways_left_even\n");
		return return_error_val();
	}
	if(k<0){
		fprintf(stderr,"ERROR k is less than 0 in "
                       "update_conways_left_even\n");
		return return_error_val();
	}
	if(get_rows(orig)!=get_rows(temp)){
		fprintf(stderr,"ERROR rows in orig is not equal to "
                       "rows in temp in "
                       "update_conways_left_even\n");
		return return_error_val();
	}
	if(get_rows(new_mat)!=get_rows(temp)){
		fprintf(stderr,"ERROR rows in new_mat is not equal to"
                       " rows in temp in "
                       "update_conways_left_even\n");
		return return_error_val();
	}
	if(get_cols(orig)!=get_cols(temp)){
		fprintf(stderr,"ERROR cols in orig is not equal to "
                       "cols in temp in "
                       "update_conways_left_even\n");
		return return_error_val();
	}
	if(get_cols(new_mat)!=get_cols(temp)){
		fprintf(stderr,"ERROR cols in new_mat is not equal "
                       "to columns in temp in "
                       "update_conways_left_even\n");
		return return_error_val();
	}
	#endif

	int row_end;
	int row_start;
	int col_end;
	int col_start;

	int ghost_left  = get_cols_left_ghost_matrix(orig);
	int ghost_right = get_cols_right_ghost_matrix(orig);
	int ghost_above = get_rows_above_ghost_matrix(orig);
	int ghost_below = get_rows_below_ghost_matrix(orig);

	if(ghost_above==0){
		row_start = 0;
	}else{
		row_start = get_row_start(orig)+1+k;
	}
	if(ghost_below==0){
		row_end = get_rows(orig);
	}else{
		row_end = get_row_end(orig)-1-k;
	}
	if(ghost_left==0){
		col_start=0;
	}else{
		col_start = get_col_start(orig)+1+k;
	}
	if(ghost_right==0){
		col_end = get_cols(orig);
	}else{
		col_end = get_col_end(orig)-1-k;
	}
	set_elems_conway_rectangle_rotation(col_start, col_end/2,
			                            row_start, row_end,
			                            orig , &temp, &new_mat, k);
	return 0;
}

int set_elems_conway_rectangle(int col_start, int col_end,
                               int row_start, int row_end,
                               const_matrix mat1, matrix * mat2,
                               int ghost){

	#ifdef _ERROR_CHECKING_ON_
	if(mat1==NULL){
		fprintf(stderr,"ERROR mat1 is NULL in "
                       "set_elems_conway_rectangle\n");
		return return_error_val();
	}
	if(mat2==NULL){
		fprintf(stderr,"ERROR mat2 is NULL in "
                       "set_elems_conway_rectangle\n");
		return return_error_val();
	}
	if(col_start<0){
		fprintf(stderr,"ERROR col_start < 0 in "
                       "set_elems_conway_rectangle\n");
		return return_error_val();
	}
	if(col_start>get_cols(mat1)){
		fprintf(stderr,"ERROR in col_start is greater "
                       "the number of cols in mat1 "
                       "set_elems_conway_rectangle\n");
		return return_error_val();
	}
	if(col_start>col_end){
		fprintf(stderr,"ERROR col_start is greater "
                       "than col_end in "
                       "set_elems_conway_rectangle\n");
		return return_error_val();
	}
	if(col_end>get_cols(mat1)){
		fprintf(stderr,"ERROR col_end is greater than"
                       " the number of cols in mat1 "
                       "in "
                       "set_elems_conway_rectangle\n");
		return return_error_val();
	}
	if(row_start<0){
		fprintf(stderr,"ERROR row_start is less than"
                       " 0 in "
                       "set_elems_conway_rectangle\n");
		return return_error_val();
	}
	if(row_start>get_rows(mat1)){
		fprintf(stderr,"ERROR row_start is greater "
                       "than the number of rows in "
                       "mat1 in "
                       "set_elems_conway_rectangle\n");
		return return_error_val();
	}
	if(row_start>row_end){
		fprintf(stderr,"ERROR row_start is greater "
                       "than row_end in "
                       "set_elems_conway_rectangle\n");
		return return_error_val();
	}
	if(row_end>get_rows(mat1)){
		fprintf(stderr,"ERROR row_end is greater than"
                       " the number of rows in mat1 "
                       "in "
                       "set_elems_conway_rectangle\n");
		return return_error_val();
	}
	if(get_rows(mat1)!=get_rows(*mat2)){
		fprintf(stderr,"ERROR number of rows in mat1"
                       " is not equal to the number "
                       "of rows in mat2 in "
                       "set_elems_conway_rectangle\n");
		return return_error_val();
	}
	if(get_cols(mat1)!=get_cols(*mat2)){
		fprintf(stderr,"ERROR number of cols in mat1"
                       "is not equal to the number "
                       "of rows in mat2 in "
                       "set_elems_conway_rectangle\n");
		return return_error_val();
	}
	#endif

	float val;
	int i;
	int j;
	for(j=col_start;j<col_end;j++){
		for(i=row_start;i<row_end;i++){
			if(init_test==0){
				conways_condition(mat1,mat2,i,j,ghost);
			}else if(init_test==1){
				val = add_all_surrounding_elems(mat1,i,j);
				set_elem(*mat2,i,j,val);
			}else if(init_test==2){
				val = get_elem(mat1,i,j)+1.0;
				set_elem(*mat2,i,j,val);
			}
		}
	}

	return 0;
}

int set_elems_conway_rectangle_rotation(int col_start, int col_end,
                                        int row_start, int row_end,
                                        const_matrix mat1 , matrix * mat2,
                                        matrix * new_mat, int k){
	
	#ifdef _ERROR_CHECKING_ON_
	if(mat1==NULL){
		fprintf(stderr,"ERROR mat1 is NULL in "
                       "set_elems_conway_rectangle_rotation\n");
		return return_error_val();
	}
	if(mat2==NULL){
		fprintf(stderr,"ERROR mat2 is NULL in "
                       "set_elems_conway_rectangle_rotation\n");
		return return_error_val();
	}
	if(col_start<0){
		fprintf(stderr,"ERROR col_start < 0 in "
                       "set_elems_conway_rectangle_rotation\n");
		return return_error_val();
	}
	if(col_start>get_cols(mat1)){
		fprintf(stderr,"ERROR in col_start is greater "
                       "the number of cols in mat1 "
                       "set_elems_conway_rectangle_rotation\n");
		return return_error_val();
	}
	if(col_start>col_end){
		fprintf(stderr,"ERROR col_start is greater "
                       "than col_end in "
                       "set_elems_conway_rectangle_rotation\n");
		return return_error_val();
	}
	if(col_end>get_cols(mat1)){
		fprintf(stderr,"ERROR col_end is greater than"
                       " the number of cols in mat1 "
                       "in "
                       "set_elems_conway_rectangle_rotation\n");
		return return_error_val();
	}
	if(row_start<0){
		fprintf(stderr,"ERROR row_start is less than"
                       " 0 in "
                       "set_elems_conway_rectangle_rotation\n");
		return return_error_val();
	}
	if(row_start>get_rows(mat1)){
		fprintf(stderr,"ERROR row_start is greater "
                       "than the number of rows in "
                       "mat1 in "
                       "set_elems_conway_rectangle_rotation\n");
		return return_error_val();
	}
		if(row_start>row_end){
		fprintf(stderr,"ERROR row_start is greater "
                       "than row_end in "
                       "set_elems_conway_rectangle_rotation\n");
		return return_error_val();
	}
	if(row_end>get_rows(mat1)){
		fprintf(stderr,"ERROR row_end is greater than"
                       " the number of rows in mat1 "
                       "in "
                       "set_elems_conway_rectangle_rotation\n");
		return return_error_val();
	}
	if(get_rows(mat1)!=get_rows(*mat2)){
		fprintf(stderr,"ERROR number of rows in mat1"
                       " is not equal to the number "
                       "of rows in mat2 in "
                       "set_elems_conway_rectangle_rotation\n");
		return return_error_val();
	}
	if(get_cols(mat1)!=get_cols(*mat2)){
		fprintf(stderr,"ERROR number of cols in mat1"
                       "is not equal to the number "
                       "of rows in mat2 in "
                       "set_elems_conway_rectangle_rotation\n");
		return return_error_val();
	}
	if(get_cols(*new_mat)!=get_cols(*mat2)){
		fprintf(stderr,"ERROR number of cols in new_mat"
                       "is not equal to the number "
                       "of rows in mat2 in "
                       "set_elems_conway_rectangle_rotation\n");
		return return_error_val();
	}
	if(get_rows(*new_mat)!=get_rows(*mat2)){
		fprintf(stderr,"ERROR number of rows in new_mat"
                       " is not equal to the number "
                       "of rows in mat2 in "
                       "set_elems_conway_rectangle_rotation\n");
		return return_error_val();
	}
	if(k<0){
		fprintf(stderr,"ERROR k is less than 0 in "
                       "set_elems_conway_rectangle_rotation\n");
		return return_error_val();
	}
	#endif

	float val;
	int i;
	int j;
	for(j=col_start;j<col_end;j++){
		for(i=row_start;i<row_end;i++){
			if(k==0 && (ghost_count)!=1){
				if(init_test==0){
					conways_condition(mat1,mat2,i,j,0);
				}else if(init_test==1){
					val = add_all_surrounding_elems(mat1,i,j);
					set_elem(*mat2,i,j,val);
				}else if(init_test==2){
					val = get_elem(mat1,i,j)+1.0;
					set_elem(*mat2,i,j,val);	
				}
			}else if(k==(ghost_count-1)){
				if(init_test==0){
					conways_condition(mat1,new_mat,i,j,0);
				}else if(init_test==1){
					val = add_all_surrounding_elems(mat1,i,j);
					set_elem(*new_mat,i,j,val);
				}else if(init_test==2){
					val = get_elem(mat1,i,j)+1.0;
					set_elem(*new_mat,i,j,val);
				}
			}else{
				if(init_test==0){
					conways_condition(mat1,mat2,i,j,0);
				}else if(init_test==1){
					val = add_all_surrounding_elems(mat1,i,j);
					set_elem(*mat2,i,j,val);
				}else if(init_test==2){
					val = get_elem(mat1,i,j)+1.0;
					set_elem(*mat2,i,j,val);	
				}
			}
		}
	}
	return 0;
}

bool conways_game(matrix * orig_mat, matrix * temp, matrix * new_mat){

	#ifdef _ERROR_CHECKING_ON_
	if(orig_mat==NULL){
		fprintf(stderr,"ERROR orig_mat is NULL in conways_game\n");
		return_error_val();
		return false;
	}
	if(new_mat==NULL){
		fprintf(stderr,"ERROR new_mat is NULL in conways_game\n");
		return_error_val();
		return false;
	}
	if((*orig_mat)->rows!=(*new_mat)->rows){
		fprintf(stderr,"ERROR orig_mat and new_mat do not have the same"
                       " number of rows in conways_game\n");
		return_error_val();
		return false;
	}
	if((*orig_mat)->cols!=(*new_mat)->cols){
		fprintf(stderr,"ERROR orig_mat and new_mat do not have the same"
                       " number of columns in conways_game\n");
		return_error_val();
		return false;
	}
	if((*orig_mat)->row_start!=(*new_mat)->row_start){
		fprintf(stderr,"ERROR orig_mat and new_mat do not have the same"
                       " number of top ghost rows in conways_game\n");
		return_error_val();
		return false;
	}
	if((*orig_mat)->row_end!=(*new_mat)->row_end){
		fprintf(stderr,"ERROR orig_mat and new_mat do not have the same"
                       " number bottom ghost rows in conways_game\n");
		return_error_val();
		return false;
	}
	if(block_on){
		if((*new_mat)->col_start!=0){
			fprintf(stderr,"ERROR new_mat has left ghost columns which "
                           "are not needed in the panel implementation "
                           "(block) in conways_game\n");
			return_error_val();
			return false;
		}
		if((*new_mat)->col_end!=(*new_mat)->cols){
			fprintf(stderr,"ERROR new_mat has right ghost columns which "
                           "are not needed in the panel implementation "
                           "(block) in conways_game\n");
			return_error_val();
			return false;
		}	
		if((*orig_mat)->col_start!=0){
			fprintf(stderr,"ERROR orig_mat has left ghost columns which "
                           "are not needed in the panel implementation "
                           "(block) in conways_game\n");
			return_error_val();
			return false;
		}
		if((*orig_mat)->col_end!=(*orig_mat)->cols){
			fprintf(stderr,"ERROR orig_mat has right ghost columns which are"
					" not needed in the panel implementation (block) "
					"in conways_game\n");
			return_error_val();
			return false;
		}
	}
	#endif

	int k;
	int temp1_flag;
	int temp2_flag;
	int new_temp_flag;
	int place_holder_flag;
	
	matrix temp1;
	matrix temp2;
	matrix new_temp;
	matrix place_holder;
	
	temp1 = *orig_mat;
	temp2 = *new_mat;
	temp1_flag = 1;
	temp2_flag = 3;

	/* Start by calculating the borders */	
	for(k=0;k<ghost_count;k++){
		if(verbose){
			printf("Updating Border ghost %d\n",k);
			fflush(NULL);
		}
		update_border_conway(temp1 , &temp2, k);
		
		if(k%2==0){
			temp1_flag = 3;
			temp2_flag = 2;
			temp1 = *new_mat;
			temp2 = *temp;
		}else{
			temp1_flag = 2;
			temp2_flag = 3;
			temp1 = *temp;
			temp2 = *new_mat;
		}
	}

	/* temp1 will be the most up to date */
	new_temp = temp1;
	new_temp_flag = temp1_flag;
    /* MPI send and receive for updating ghost
     * rows can be made */
	#ifdef _MPI_H_
	if(block_on==false){
		/* For the row block config we don't need
           to send the columns at all */
		if(!(my_rank==(num_proc-1) && my_rank==0)){
			send_recv_ghost_cols(&new_temp);
		}
		/* We can not asynchronously send the rows 
		 * at this point if using the checker config
		 * this is because the col send recv has to
		 * be completed before the rows can be updated
		 */
	}
	#endif

	/* This means that temp1 currenty has the most up to 
     * date data, new_temp points to the same matrix we
     * want to resuse temp1 but without affecting the 
     * most upto date matrix so we point temp1 to a
     * different matrix. */
	if(temp1_flag==2){
		
		temp1_flag = 1;
		temp1      = *orig_mat;
		temp2_flag = 3;
		temp2      = *new_mat;	
	}else if(temp1_flag==3){
		temp1_flag = 1;
		temp1      = *orig_mat;
		temp2_flag = 2;
		temp2      = *temp;
	}

	if(verbose){
		if(init_test==0){
			if(block_on){
				printf("After sending and receiving ghost rows\n");
			}else{
				printf("After sending and receiving ghost cols\n");
			}
			print_matrix(new_temp);
			fflush(NULL);
		}else{
			if(block_on){
				printf("After sending and receiving ghost rows\n");
			}else{
				printf("After sending and receiving ghost cols\n");
			}
			print_matrix_compressed(new_temp);
			fflush(NULL);
		}
	}

	/* Now the inner ghost rows and columns can be updated */
	/* Split the matrix up into two parts this is so we can 
     * take advantage of asynchrounous communication. There
     * exists two different scenarios the case where the 
     * number of ghost rows/columns is odd or where they are 
     * even. If they are odd then we will split the matrix 
     * down the middle. */
	for(k=0;k<ghost_count;k++){

		update_conways_left(temp1,temp2, new_temp, k);

		if(k==ghost_count/2 && ghost_count%2==1){
			/* This means that we are in the middle of the 
			 * calculation. if it is an odd number of ghost cols
			 * and rows */
			#ifdef _MPI_H_
			if(!(my_rank==(num_proc-1)&&my_rank==0)){
				if(send_recv_sync==false){
					/* Halfway through the rest of the matrix is the optimal
					 * time to send the rows if using the checker scheme. 
					 * however we first have to make sure the column message
					 * passing is complete
					 */
					if(block_on==false){
						if(my_grid_col==0){
							MPI_Wait(&send_request3,&status);
						}else if (my_grid_col==grid_dim_col-1){
							MPI_Wait(&recv_request3,&status);
						}else{
							MPI_Wait(&recv_request3,&status);
							MPI_Wait(&send_request3,&status);
						}
						if(my_grid_col==0){
							MPI_Wait(&recv_request4,&status);
						}else if (my_grid_col==grid_dim_col-1){
							MPI_Wait(&send_request4,&status);
						}else{
							MPI_Wait(&recv_request4,&status);
							MPI_Wait(&send_request4,&status);
						}
					}
				} 
				/* Regardless of wehther this is row block or
				 * checker we need to send the rows */
				send_recv_ghost_rows(&new_temp);
				if(verbose){
					if(init_test==0){
						if(block_on==false){
							printf("After sending and receiving ghost rows"
									" if checkboard\n");
							printf("Should be halfway through core matrix "
									"calculation\n");
							printf("Here is the matrix that should have "
									"received the rows if synchronous and odd\n");
							print_matrix(new_temp);
							printf("Here is the matrix with half the core"
									" calculated\n");
							printf("Temp2 flag %d\n",temp2_flag);
							print_matrix(temp2);
							fflush(NULL);
						}
					}else{
						if(block_on==false){
							printf("After sending and receiving ghost rows"
									" if checkerboard \n");
							printf("Should be halfway through core matrix "
									"calculation\n");
							printf("Here is the matrix that should have "
									"received the rows if synchronous and odd\n");
							print_matrix_compressed(new_temp);
							printf("Here is the matrix with half the core"
									" calculated\n");
							printf("Temp2 flag %d\n",temp2_flag);
							print_matrix_compressed(temp2);
							fflush(NULL);
						}
					}
				}



			}
			#endif
		}

		/* For the right side of the matrix */
		update_conways_right(temp1,temp2,new_temp,k);
		place_holder_flag = temp1_flag;
		place_holder      = temp1;
		temp1_flag        = temp2_flag;
		temp1             = temp2;
		temp2_flag        = place_holder_flag;
		temp2             = place_holder;

		if(k==(ghost_count/2-1) && ghost_count%2==0){
			/* This means that we are in the middle of the 
			 * calculation. if it is an even number of ghost cols
			 * and rows */
			#ifdef _MPI_H_
			if(!(my_rank==(num_proc-1)&&my_rank==0)){
				if(send_recv_sync==false){
					/* Halfway through the rest of the matrix is the optimal
					 * time to send the rows if using the checker scheme. 
					 * however we first have to make sure the column message
					 * passing is complete
					 */
					if(block_on==false){
						if(my_grid_col==0){
							MPI_Wait(&send_request3,&status);
						}else if (my_grid_col==grid_dim_col-1){
							MPI_Wait(&recv_request3,&status);
						}else{
							MPI_Wait(&recv_request3,&status);
							MPI_Wait(&send_request3,&status);
						}
						if(my_grid_col==0){
							MPI_Wait(&recv_request4,&status);
						}else if (my_grid_col==grid_dim_col-1){
							MPI_Wait(&send_request4,&status);
						}else{
							MPI_Wait(&recv_request4,&status);
							MPI_Wait(&send_request4,&status);
						}
					}
				} 
				/***************CHANGING FROM temp1 ***********/
				send_recv_ghost_rows(&new_temp);
				if(verbose){
					if(init_test==0){
						if(block_on==false){
							printf("After sending and receiving ghost rows"
									" if checkboard\n");
							printf("Should be halfway through core matrix "
									"calculation\n");
							printf("Here is the matrix that should have "
									"received the rows if even and syhncronous\n");
							print_matrix(new_temp);
							printf("Here is the matrix with half the core"
									" calculated\n");
							print_matrix(temp1);
							fflush(NULL);
						}
					}else{
						if(block_on==false){
							printf("After sending and receiving ghost rows"
									" if checkerboard \n");
							printf("Should be halfway through core matrix "
									"calculation\n");
							printf("Here is the matrix that should have "
									"received the rows if even and synchronous\n");
							print_matrix_compressed(new_temp);
							printf("Here is the matrix with half the core"
									" calculated\n");
							print_matrix_compressed(temp1);
							fflush(NULL);
						}
					}
				}

			}
			#endif
		}
	}
	#ifdef _MPI_H_
	if(!(my_rank==(num_proc-1)&&my_rank==0)){

		if(send_recv_sync==false){
		
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
		}
	}
	if(verbose){
		printf("Printing new_temp after completion of"
               "conway's game\n");
		fflush(NULL);
		if(init_test==0){
			print_matrix(new_temp);
		}else{
			print_matrix_compressed(new_temp);
		}
	}
	#endif

	/* Update which matrix is the newest matrix */
	*new_mat = new_temp;

	/* And the oldest matrix */
	*orig_mat = temp2;

	/* Middle matrix */
	*temp     = temp1;	

	if(verbose){
		printf("After conways completed\n");
		if(init_test==0){
			print_matrix(new_temp);
			fflush(NULL);
		}else{
			print_matrix_compressed(new_temp);
			fflush(NULL);
		}
	}

	return 0;
}

float add_all_surrounding_elems(const_matrix orig_mat, int i, int j){

	#ifdef _ERROR_CHECKING_ON_
	if(orig_mat==NULL){
		fprintf(stderr,"ERROR *orig_mat is NULL in "
                       "add_all_surrounding_elems\n");
		return (float) return_error_val();
	}
    if(i<0){
        fprintf(stderr,"ERROR row is less than 0 in "
                       "add_all_surrounding_elems\n");
		return (float) return_error_val();
    }
    if(j<0){
        fprintf(stderr,"ERROR col is less than 0 in "
                       "add_all_surrounding_elems\n");
		return (float) return_error_val();
    }
    if(i>=orig_mat->rows){
        fprintf(stderr,"ERROR row greater than rows in matrix "
                       "add_all_surrounding_elems\n");
		return (float) return_error_val();
    }
    if(j>=orig_mat->cols){
        fprintf(stderr,"ERROR col greater than cols in matrix "
                       "add_all_surrounding_elems\n");
		return (float) return_error_val();
    }
	#endif
	float sum = 0.0;
    float value;
    float NW;
    float N;
    float NE;
    float E;
    float W;
    float SW;
    float S;
    float SE;

    value = get_elem(orig_mat,i,j);
    if(within_matrix(orig_mat,i-1,j-1)){
        NW = get_elem(orig_mat,i-1,j-1);
        sum += NW;
    }
    if(within_matrix(orig_mat,i,j-1)){
        N = get_elem(orig_mat,i,j-1);
        sum += N;
    }
    if(within_matrix(orig_mat,i+1,j-1)){
        NE = get_elem(orig_mat,i+1,j-1);
        sum += NE;
    }
    if(within_matrix(orig_mat,i-1,j)){
        W = get_elem(orig_mat,i-1,j);
        sum += W;
    }
    if(within_matrix(orig_mat,i+1,j)){
        E = get_elem(orig_mat,i+1,j);
        sum += E;
    }
    if(within_matrix(orig_mat,i-1,j+1)){
        SW = get_elem(orig_mat,i-1,j+1);
        sum += SW;
    }
    if(within_matrix(orig_mat,i,j+1)){
        S = get_elem(orig_mat,i,j+1);
        sum += S;
    }
	if(within_matrix(orig_mat,i+1,j+1)){
		SE = get_elem(orig_mat,i+1,j+1);
		sum += SE;
	}

	return sum;

}
