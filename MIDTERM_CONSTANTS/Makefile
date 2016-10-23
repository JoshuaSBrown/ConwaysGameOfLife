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
CFAGS = -Wall -Wextra -std=c99 -pedantic -g --save-temps

ALL : test_midterm_constants 

test_midterm_constants : test_midterm_constants.c midterm_constants.o
	$(CC) $(CFLAGS) -o test_midterm_constants test_midterm_constants.c midterm_constants.o -lm 

midterm_constants.o : ../MIDTERM_CONSTANTS/midterm_constants.c 
	$(CC) $(CFLAGS) -c ../MIDTERM_CONSTANTS/midterm_constants.c

.PHONY : clean
clean :
	$(RM) *.o
