/**************************************/
/*          Midterm 3
 *          Joshua Brown
 *          Scientific Computing
 *          midterm_constants.c function
 *          Conway's Game of Life
 */
/**************************************/

#ifndef _MIDTERM_IO_H_
#define _MIDTERM_IO_H_

#include <stdbool.h>

#include "../MIDTERM_MATRIX/midterm_matrix.h"

/* This function looks in a string for:
 *  1 a single interger
 *  2 a list of intergers 
 *  3 a range of intergers
 *
 * The format of these 3 types is as shown
 * below:
 * Single "5" 
 * List   "434,692, 12345"
 * range  "23-235"
 *
 * Warning lists and ranges can not be combined
 * you cannot have 
 * "234,1-54"
 *
 * You can not have more than one dash as well
 * such as 
 * "245-963-3245"
 * 
 * The return value is an interger:
 *  o -1 for malformed input or invalid string layout
 *  o  0 Single value input
 *  o  1 List   comma delliminated input
 *  o  2 Range  dash deliminated input
 */
int single_list_range(char * str);

/* This function will check for the number of dashs
 * "-" that are within a string. The return value is
 *  o -1 Malformed string input
 *  o  0 No dashs
 *  o >0 A value greater than 0 indicating the number of
 *       dashes
 */
int string_check_for_dashs(char * str);

/* Searches the string for commas "," returns:
 *  o -1 Malformed input
 *  o  0 No commas
 *  o >0 A value greater than 0 indicating the number of
 *       commas
 */
int string_check_for_commas(char * str);

/* This function basically checks that the number of dashs
 * and commas are appropriate and that the string follows 
 * the correct rules. Returns:
 *  o -1 If malformed input or string in incorrectly 
 *       formatted.
 *  o  0 If there are no commas or dashes indicated there
 *       is only a single interger value.
 *  o >0 A value indicating the number of commas or dashs
 */
int pgm_input_number_values(char *str);

/* Takes a string that is correctly formatted with either
 * commas or a dash and reads the interger values into a 
 * matrix. Returns
 *  o -1 If malformed input string or matrix
 *  o  0 If successful 
 *
 * Warning the matrix must already be malloced before it
 * is passed to this function */
int pgm_read_values_to_array(matrix mat, char * str);

/* This sorts the command line options */
int options(int argc, char **argv);

/* Prints a message to the screen indicating how to use the
 * program */
void print_help_message(void);

/* Reads the pgm file to matrix mat_a and mat_b */
bool read_pgm( matrix * mat_a, matrix * mat_b);

/* Writes a matrix to a .pgm file count is an interger used 
 * to name the .pgm file. */
bool write_pgm( const_matrix mat_a, int count);

/* Defines a matrix based on the constants defined in 
 * midpoint_constants.h */
matrix define_matrix_blockchecker(void);

double standard_deviation(double mean, int iteration);

int calculate_statistics(int iteration);
#endif
