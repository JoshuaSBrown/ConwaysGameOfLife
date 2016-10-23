#include <stdio.h>
#include <float.h>

#include "mpi.h"

int main(int argc, char *argv[]){


	MPI_File fh;
	int buf[10], rank;
	char buf2[10];
	char buf3[] = "P5\n900 900\n255\n"; 
	int siz = sizeof(buf3);
	int i;
	for(i=0;i<10;i++){
		buf[i]=255;
		buf2[i]=255;
	}
	MPI_Init(0,0);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_File_open(MPI_COMM_WORLD,"test.out",
	              MPI_MODE_CREATE | MPI_MODE_WRONLY,
                  MPI_INFO_NULL, &fh);
	if(rank==0){
		MPI_File_write(fh,buf3,siz,MPI_CHAR,MPI_STATUS_IGNORE);
		MPI_File_write(fh,buf2,10,MPI_CHAR,MPI_STATUS_IGNORE);
	}
	MPI_File_close(&fh);

	float buf4[1];
	buf4[0] = FLT_MAX;
	MPI_File fh2;
	MPI_File_open(MPI_COMM_WORLD,"testa.out",
                  MPI_MODE_CREATE | MPI_MODE_WRONLY,
                  MPI_INFO_NULL, &fh2);

	if(rank==0){
		MPI_File_write(fh2,buf4,1,MPI_FLOAT,MPI_STATUS_IGNORE);
		MPI_File_write(fh2,buf4,1,MPI_CHAR,MPI_STATUS_IGNORE);
	}	

	MPI_File_close(&fh2);

	float buf5[1];
	buf5[0] = FLT_MIN;
	MPI_File fh3;
	MPI_File_open(MPI_COMM_WORLD,"testb.out",
                  MPI_MODE_CREATE | MPI_MODE_WRONLY,
                  MPI_INFO_NULL, &fh3);

	if(rank==0){
		MPI_File_write(fh3,buf5,1,MPI_FLOAT,MPI_STATUS_IGNORE);
	}	

	MPI_File_close(&fh3);


	MPI_Finalize();
	return 0;
}
