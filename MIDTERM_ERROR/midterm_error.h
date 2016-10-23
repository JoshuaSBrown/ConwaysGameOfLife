/**************************************/
/*          Midterm 3
 *          Joshua Brown
 *          Scientific Computing
 *          midterm_error.h function
 *          Conway's Game of Life
 */
/**************************************/
#ifndef _MIDTERM_ERROR_H_
#define _MIDTERM_ERROR_H_

#ifndef _MPI_H_
#define _MPI_H_
#endif
#ifndef _ERROR_CHECKING_ON_
#define _ERROR_CHECKING_ON_
#endif
#ifndef _FORCE_HARD_CRASH_
#define _FORCE_HARD_CRASH_
#endif

/* Used to write stdout to a file
*/
static int fd;

/* Used to keep track of the position
 * in the file we are writing too. */
static fpos_t pos;

/* Filename used to output standard out */
char filename_out[100];

/* Initialized file Name */
void create_filename_out(int num);

/* The following two functions were taken
 * from http://stackoverflow.com/questions/14543443/in-c-how-do-you-redirect-stdin-stdout-stderr-to-files-when-making-an-execvp-or
 */
void switch_stdout(const char *newStream);

void revert_stdout();

/* Prints ERROR message indicating that
 * the pointer described with ptr_name 
 * in the function described in 
 * function_name is NULL. The length
 * variable describeds how long the
 * combined two characters arrays are
 */
void print_NULL_error(char * ptr_name,
    char * function_name, int length);

/* Prints ERROR message indicating that
 * the interger value described by 
 * var_name is negative when it should
 * be positive in the function
 * function_name. The length descirbes
 * how long the combined character 
 * arrays are.
 */
void print_negative_int_error(char * var_name,
    char * function_name, int length);

/* Function will exit correctly depending
 * on whether it is being run in parallel
 * or not. If it is not forced to exit will
 * return a value of -1. 
 */
int return_error_val();

#endif
