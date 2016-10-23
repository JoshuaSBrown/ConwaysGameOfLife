# /**************************************/
# /*          Midterm 2
#  *          Joshua Brown
#  *          Scientific Computing
#  *          Makefile for test_midterm_matrix.c
#  *          Conway's Game of Life
#  */
# /**************************************/

CC = mpicc
CFAGS = -Wall -Wextra -std=c99 -pedantic -g --save-temps

ALL : test_midterm_io create_pgm

test_midterm_io : test_midterm_io.c midterm_io.o midterm_error.o midterm_pprintf.o midterm_matrix.o midterm_constants.o
	$(CC) $(CFLAGS) -o test_midterm_io test_midterm_io.c midterm_io.o midterm_matrix.o midterm_error.o  midterm_pprintf.o midterm_constants.o -lpopt -lm

create_pgm : create_pgm.c midterm_io.o midterm_constants.o midterm_matrix.o midterm_pprintf.o midterm_error.o
	$(CC) $(CFLAGS) -o create_pgm create_pgm.c midterm_io.o midterm_constants.o midterm_matrix.o midterm_pprintf.o midterm_error.o -lpopt -lm

midterm_io.o : midterm_io.c 
	$(CC) $(CFLAGS) -c midterm_io.c -lpopt

midterm_pprintf.o : ../MIDTERM_PPRINTF/midterm_pprintf.c
	$(CC) $(CFLAGS) -c ../MIDTERM_PPRINTF/midterm_pprintf.c

midterm_matrix.o : ../MIDTERM_MATRIX/midterm_matrix.c 
	$(CC) $(CFLAGS) -c ../MIDTERM_MATRIX/midterm_matrix.c -lm

midterm_error.o : ../MIDTERM_ERROR/midterm_error.c  
	$(CC) $(CFLAGS) -c ../MIDTERM_ERROR/midterm_error.c 

midterm_constants.o : ../MIDTERM_CONSTANTS/midterm_constants.c
	$(CC) $(CFLAGS) -c ../MIDTERM_CONSTANTS/midterm_constants.c -lm

.PHONY : clean
clean :
	$(RM) *.o
