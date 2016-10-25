#!/bin/bash
#SBATCH --output=out2.txt
#SBATCH --partition=development
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=9
#SBATCH -t 00:02:00
#SBATCH --mail-user=joshbro42867@yahoo.com

cd $WORK
rm -rf *
mkdir midterm3
cd midterm3/
cp $HOME/ConwaysGameOfLife/MIDTERM_MAIN/midterm_main .
cp $HOME/ConwaysGameOfLife/MIDTERM_MAIN/conways_input.pgm .
ibrun -np 9 ./midterm_main -f 2700x2700.pgm -i 4 -g -s > midterm3_out.txt
cp *txt $HOME/ConwaysGameOfLife/MIDTERM_MAIN/
cp *out $HOME/ConwaysGameOfLife/MIDTERM_MAIN/
cd $HOME
