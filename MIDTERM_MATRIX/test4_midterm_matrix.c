#include <stdio.h>

#include "midterm_matrix.h"

int main(void){

	char file_name[] = "matlab_data_file.txt";

	matrix mat = new_matrix(10,10);
	set_elem(mat,4,4,1.0);
	set_elem(mat,3,3,1.0);
	set_elem(mat,3,4,1.0);
	set_elem(mat,5,4,1.0);
	print_matrix(mat);
	generate_data_file(file_name,mat);
		
	delete_matrix(&mat);
	return 0;
}
