#include "stdio.h"
#include "midterm_io.h"
#include "../MIDTERM_CONSTANTS/midterm_constants.h"
#include "mpi.h"

int main(){

	MPI_Init(0,0);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

	initialize_grid_checker();

	int dim = 6;
	sub_matrix_dim_columns = dim*2;
	sub_matrix_dim_rows    = dim*2;

	matrix mat = new_matrix(sub_matrix_dim_rows,
                            sub_matrix_dim_columns);
	
//	set_elem(mat,3,3,1.0);
//	set_elem(mat,3,4,1.0);
//	set_elem(mat,3,2,1.0);

	int i;
	int j;

	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++){
			set_elem(mat,i,j,1.0);
		}
	}

	for(i=dim;i<dim*2;i++){
		for(j=dim;j<dim*2;j++){
			set_elem(mat,i,j,1.0);
		}
	}
	print_matrix(mat);
	write_pgm( mat, 99);

	MPI_Finalize();
}
