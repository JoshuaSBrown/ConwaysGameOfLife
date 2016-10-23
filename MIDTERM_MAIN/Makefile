# /**************************************/
# /*          Midterm 3
#  *          Joshua Brown
#  *          Scientific Computing
#  *          Makefile for midterm_matrix.c 
#  *               and for test_midterm_matrix.c
#  *          Conway's Game of Life
#  */
# /**************************************/

CC = mpicc
CFAGS = -Wall -Wextra -std=c99 -O3 -pedantic -g -static

ALL : midterm_main test_midterm_main 

midterm_main : midterm_main.c midterm_matrix.o midterm_error.o midterm_constants.o midterm_io.o midterm_pprintf.o
	$(CC) $(CFLAGS) -o midterm_main midterm_main.c midterm_matrix.o midterm_error.o  midterm_constants.o midterm_io.o midterm_pprintf.o -lpopt -lm 

test_midterm_main : test_midterm_main.c midterm_matrix.o midterm_error.o midterm_constants.o midterm_io.o midterm_pprintf.o
	$(CC) $(CFLAGS) -o test_midterm_main test_midterm_main.c midterm_matrix.o midterm_error.o  midterm_constants.o midterm_io.o midterm_pprintf.o -lpopt -lm 

midterm_matrix.o : ../MIDTERM_MATRIX/midterm_matrix.c midterm_error.o 
	$(CC) $(CFLAGS) -c ../MIDTERM_MATRIX/midterm_matrix.c midterm_error.o 

midterm_constants.o : ../MIDTERM_CONSTANTS/midterm_constants.c 
	$(CC) $(CFLAGS) -c ../MIDTERM_CONSTANTS/midterm_constants.c

midterm_io.o : ../MIDTERM_IO/midterm_io.c
	$(CC) $(CFLAGS) -c ../MIDTERM_IO/midterm_io.c -lm

midterm_pprintf.o : ../MIDTERM_PPRINTF/midterm_pprintf.c
	$(CC) $(CFLAGS) -c ../MIDTERM_PPRINTF/midterm_pprintf.c

midterm_error.o : ../MIDTERM_ERROR/midterm_error.c  
	$(CC) $(CFLAGS) -c ../MIDTERM_ERROR/midterm_error.c 

.PHONY : clean
clean :
	$(RM) *.o
