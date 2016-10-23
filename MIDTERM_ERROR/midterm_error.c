/**************************************/
/*          Midtermp 3 
 *          Joshua Brown
 *          Scientific Computing
 *          midterm_error.c function
 *          Conway's Game of Life
 */
/**************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "midterm_error.h"

#ifdef _MPI_H_
#include "mpi.h"
#endif

void create_filename_out(int num){
	sprintf(filename_out,"file_num%d.out",num);
}

void switch_stdout(const char *newStream){
	fflush(stdout);
	fgetpos(stdout,&pos);
	fd = dup(fileno(stdout));
	freopen(newStream,"w",stdout);
}

void revert_stdout(){
	fflush(stdout);
	dup2(fd, fileno(stdout));
	close(fd);
	clearerr(stdout);
	fsetpos(stdout, &pos);
}

void print_NULL_error(char * ptr_name, char * function_name, int length){
    int num_char = length + 22;
    char message[num_char];
    sprintf(message,"ERROR %s is NULL in %s\n",ptr_name,function_name);
    printf("%s",message);
}

void print_negative_int_error(char * var_name, char * function_name, int length){
    int num_char = length + 26;
    char message[num_char];
    sprintf(message,"ERROR %s is negative in %s\n",var_name,function_name);
    printf("%s",message);
}

int return_error_val(){
	#ifdef _FORCE_HARD_CRASH_
	#ifdef _MPI_H_
	MPI_Abort(MPI_COMM_WORLD,-1);
	#else
	exit(1);
	#endif
	#endif
	return -1;
}
