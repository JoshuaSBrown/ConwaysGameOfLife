#include <stdio.h>

#include "mpi.h"

int main(int argc, char *argv[]){


	MPI_File fh;
	MPI_Offset offset;
	int  rank;
	MPI_Init(0,0);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	char buf2[15];	
	if(rank==0){	 
		sprintf(buf2, "P5\n900 900\n255\n"); 

	}else{
		sprintf(buf2, "P6\n500 500\n122\n"); 
	}
	int siz = sizeof(buf2);
	offset = siz*rank; 
	MPI_File_open(MPI_COMM_WORLD,"test4.out",
	              MPI_MODE_CREATE | MPI_MODE_WRONLY,
                  MPI_INFO_NULL, &fh);
	
	MPI_File_seek(fh, offset, MPI_SEEK_SET);
	MPI_File_write(fh,buf2,siz,MPI_CHAR,MPI_STATUS_IGNORE);
	MPI_File_close(&fh);
	MPI_Finalize();
	return 0;
}
