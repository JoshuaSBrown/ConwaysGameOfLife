/**************************************/
/*          Midterm 3
 *          Joshua Brown
 *          Scientific Computing
 *          midterm_io.c function
 *          Conway's Game of Life
 */
/**************************************/ 
/*
 * The pgm_read code was initially provided by Micheal Oberg and
 * has beeen edited by Joshua Brown. 
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <popt.h>
#include <stdbool.h>
#include <math.h>

#include "../MIDTERM_CONSTANTS/midterm_constants.h"
#include "../MIDTERM_ERROR/midterm_error.h"
#include "../MIDTERM_MATRIX/midterm_matrix.h"
#include "../MIDTERM_PPRINTF/midterm_pprintf.h"
#include "midterm_io.h"

#ifdef _MPI_H_
#include "mpi.h"
#endif
/* Will let the user know if string contains a single iteration or
 * if we will print out a list of values, or a range of iterations
 * Will return:
 *   -1 Malformed string input
 *    0 Single value input
 *    1 List comma deliminated
 *    2 range dash deliminated
 */
int single_list_range(char * str){

	#ifdef _ERROR_CHECKING_ON_
	if(str==NULL){
		fprintf(stderr,"ERROR str is NULL in single_list_range\n");
		return return_error_val();
	}
	#endif

    int i;
    int dash;
    int comma;
    int length;
    int count_dash;

    dash = 0;
    comma = 0;
    count_dash = 0;
    length = strlen(str);

    for(i=0;i<=length;i++){
        if(str[i]==','){
            comma=1;
            break;
        }
    }

    for(i=0;i<=length;i++){
        if(str[i]=='-'){
            count_dash++;
            dash=1;
            if(count_dash>1){
				fprintf(stderr,"ERROR only allowed one '-' in the string"
                               "when specifying a range\n");
                return -1;
            }
        }
    }

    if(dash==1 && comma==1){
		fprintf(stderr,"ERROR malformed input cannot provide a string" 
                       " with both commas and dashs. Choose one or the"
                       " other\n");
        return -1;
    }
    if(dash==1){
        return 2;
    }
    if(comma==1){
        return 1;
    }
    if(length!=0){
        /* string is not empty but no commas or dashs */
        return 0;
    }
	fprintf(stderr,"ERROR string was empty\n");
    return -1;
}

int string_check_for_dashs(char * str){

	#ifdef _ERROR_CHECKING_ON_
    if(!str){
        fprintf(stderr,"ERROR string is NULL in string_check_for_dashs\n");
        return return_error_val();
    }
	#endif

    int length = strlen(str);
    int i;
    int count_dashs = 0;
    for(i=0;i<=length;i++){
		printf("%c",str[i]);
        if(str[i]=='-'){
            count_dashs++;
        }
    }
    return count_dashs;
}

int string_check_for_commas(char * str){

	#ifdef _ERROR_CHECKING_ON_
    if(!str){
        fprintf(stderr,"ERROR str is NULL in string_check_for_commas\n");
        return return_error_val();
    }
	#endif

    int length = strlen(str);
    int i;
    int count_commas = 0;
    for(i=0;i<=length;i++){
        if(str[i]==','){
            count_commas++;
        }
    }
    return count_commas;
}


int pgm_input_number_values( char * str){

	#ifdef _ERROR_CHECKING_ON_
	if(str==NULL){
		fprintf(stderr,"ERROR str is NULL in pgm_input_number_values\n");
		return return_error_val();
	}
	#endif
    int num_commas = string_check_for_commas(str);
    int num_dashs  = string_check_for_dashs(str);

    char * p;
    int i;

    if(num_commas!=0 && num_dashs!=0){
        fprintf(stderr,"ERROR When printing pgm files either specify a range "
                       "in the following format:\n"
                       "  \"value1-value2\" \n"
                       "or\n"
                       "Specify the values of interest comma seperated:\n"
                       "  \"value1,value2,value3\"\n");
        return return_error_val();
    }

    if(num_commas==0 && num_dashs==0){
    /* This means we should only have one number */
        return 0;
    }
    if(num_commas>0){
        return num_commas;
    }
    if(num_dashs>0){
        return num_dashs;
    }

    return -1;
}


int pgm_read_values_to_array( matrix mat, char *str){

	#ifdef _ERROR_CHECKING_ON_
    if(!mat){
        fprintf(stderr,"ERROR mat is NULL in pgm_read_values\n");
        return return_error_val();
    }
    if(!str){
        fprintf(stderr,"ERROR str is NULL in pgm_read_values\n");
        return return_error_val();
    }
	#endif

    int num_commas = string_check_for_commas(str);
    int num_dashs  = string_check_for_dashs(str);
    char * token;
    char * str_copy = strdup(str);
    char * str_copy2 = strdup(str);
    int i = 0;

    printf("\nNumber of commas %d string %s\n",num_commas,str);
    if(num_commas>0){
        token = strsep( &str_copy,",");
        printf("\n%s\n",token);
        set_elem(mat,i,0,(float)atoi(token));
        for(i=1;i<=num_commas;i++){
            token = strsep( &str_copy,",");
            if(token==NULL){
                break;
            }
            printf("%s\n",token);
            set_elem(mat,i,0,(float)atoi(token));
        }
        return 0;
    }


    i = 0;
    if(num_dashs>0){
        token = strsep( &str_copy2,"-");
        printf("\n%s\n",token);
        set_elem(mat,i,0,(float)atoi(token));
        for(i=1;i<=num_dashs;i++){
            token = strsep( &str_copy2,"-");
            if(token==NULL){
                break;
            }
            set_elem(mat,i,0,(float)atoi(token));
        }

		#ifdef _ERROR_CHECKING_ON_
		if(get_elem(mat,0,0)>get_elem(mat,1,0)){
			fprintf(stderr,"ERROR When entering a range the smaller number"
                           " must come first\n");
			return return_error_val();
		}
		#endif
        return 0;
    }

    /* Single interger no dashs or commas */
    set_elem(mat,0,0,(float)atoi(str));
    free(str_copy);
    free(str_copy2);
    return 0;
}

int options(int argc, char **argv){

	/* Local Flag */
	bool flag_file = false;

	/* Initialize Values */
	verbose        = false;
	block_on       = true;
	send_recv_sync = false;
	init_data      = false;
	init_matrix_increment_pgm = NULL; 
	flag_increment_pgm = -1; 
	init_file = (char *) malloc(sizeof(char)*50);
	
	/* Initializeing matrix so that it will
     * divide up among the processors evenly
     * with a two by two matrix for each proc
     */
	init_rows   = num_proc*2;
	init_cols   = num_proc*2;
	init_iter   = 10;
	init_print  = 0;
	init_test   = 0;
	ghost_count = 1;
	poptContext POPT_Context;
	char POPT_Ret;

	int array_length;

	struct poptOption optionsTable[] =
	{
		{ "verbose"      , 'v' , POPT_ARG_NONE   , 0               , 1 , "Verbose Output"                  , 0 },
		{ "rows"         , 'r' , POPT_ARG_INT    , &init_rows      , 2 , "Number of rows of matrix"        , 0 },
		{ "cols"         , 'c' , POPT_ARG_INT    , &init_cols      , 3 , "Number of cols of matrix"        , 0 },
		{ "block"        , 'b' , POPT_ARG_NONE   , 0               , 4 , "Use block implementation"        , 0 },
		{ "checkerboard" , 'g' , POPT_ARG_NONE   , 0               , 5 , "Use checkerboard implementation" , 0 },
		{ "synchronous"  , 's' , POPT_ARG_NONE   , 0               , 6 , "Use synchronous communication"   , 0 },
		{ "asynchronous" , 'a' , POPT_ARG_NONE   , 0               , 7 , "Use asynchronous communication"  , 0 },
		{ "iterations"   , 'i' , POPT_ARG_INT    , &init_iter      , 8 , "Number of iterations"            , 0 },
		{ "print"        , 'p' , POPT_ARG_INT    , &init_print     , 9 , "Number of iterations per print " , 0 },
		{ "data"         , 'd' , POPT_ARG_NONE   , 0               , 10, "Output data file for matlab "    , 0 },
		{ "print-pgm"    , 'm' , POPT_ARG_STRING , &init_print_pgm , 11, "What iterations to print pgm   " , 0 },
		{ "ghost_count"  , 'e' , POPT_ARG_INT    , &ghost_count    , 12, "Number of ghost rows and cols"   , 0 },
		{ "test"         , 't' , POPT_ARG_INT    , &init_test      , 13, "Run test number"                 , 0 },
		{ "file"         , 'f' , POPT_ARG_STRING , &init_file      , 14, "Input file"                      , 0 },
		{ "help"         , 'h' , POPT_ARG_NONE   , 0               , 15, "Help"                            , 0 },
		POPT_AUTOHELP
		{NULL, '\0', 0, NULL, 0}
	};

	POPT_Context = poptGetContext(NULL, argc, (const char**) argv, optionsTable, 0);
	poptSetOtherOptionHelp(POPT_Context, "[ Try --help for a more detailed description of the options]");

		while((POPT_Ret = poptGetNextOpt(POPT_Context)) >=0){
			switch(POPT_Ret){
				
				case 1:
					printf("Verbose output set.\n");
					fflush(NULL);
					verbose = true;
					break;
			case 2:
				printf("Rows set %d.\n",init_rows);
				fflush(NULL);
				break;
			case 3:
				printf("Cols set %d.\n",init_cols);
				fflush(NULL);
				break;
			case 4:
				printf("Using block algorithm.\n");
				block_on = true;
				fflush(NULL);
				break;
			case 5:
				printf("Using checkerboard algorithm.\n");
				block_on = false;
				fflush(NULL);
				break;
			case 6:
				printf("Using synchronous communication.\n");
				send_recv_sync = true;
				fflush(NULL);
				break;
			case 7:
				printf("Using asynchronous communication.\n");
				fflush(NULL);
				break;
			case 8:
				printf("Setting Iterations %d.\n",init_iter);
				fflush(NULL);
				break;
			case 9:
				printf("Setting print iterations %d.\n",init_print);
				fflush(NULL);
				break;
			case 10:
				printf("Setting matlab data file print\n");
				if(num_proc>1){
					fprintf(stderr,"ERROR only setup to print data file"
							" if running\nin serial.\n");
					return_error_val();
				}
				init_data = true;
				fflush(NULL);
				break;
			case 11:
				printf("Setting print pgm on iterations %s.\n",init_print_pgm);
				flag_increment_pgm = single_list_range(init_print_pgm);
				array_length = pgm_input_number_values(init_print_pgm);
				init_matrix_increment_pgm = new_matrix(array_length+1,1);
				pgm_read_values_to_array(init_matrix_increment_pgm , init_print_pgm);
				fflush(NULL);
				break;
			case 12:
				printf("Setting ghost count %d.\n",ghost_count);
				fflush(NULL);
				if(ghost_count<0){
					fprintf(stderr,"ERROR Must have a ghost count of at"
							" least 1\n");
					return_error_val();	
				}
				break;
			case 13:
				printf("Setting test number %d.\n",init_test);
				fflush(NULL);
				break;
			case 14:
				printf("Setting File %s.\n",init_file);
				fflush(NULL);
				flag_file = true;
				break;
			case 15:
				if(my_rank==0){
					print_help_message();
					fflush(NULL);
				}
				break;
		}
	}

	if(verbose && my_rank==0){
		printf("\n");
		printf("Parameter status if values were not set will use default.\n");
		printf("Verbose option %d\n", verbose     );
		printf("Rows           %d\n", init_rows   );
		printf("Cols           %d\n", init_cols   );
		printf("Iterations     %d\n", init_iter   );
		printf("Ghost count    %d\n", ghost_count );
	}

	/* Check inputs */
	if(ghost_count>1 && init_iter%ghost_count!=0){
		fprintf(stderr,"ERROR The number of iterations does not "
				"divide evenly among the ghost_count.\n");
		return return_error_val();
	}
	if(ghost_count>1 && init_print%ghost_count){
		fprintf(stderr,"ERROR The print iterations does not divide up"
				"evenly among the ghost_count.\n");
		return return_error_val();
	}
	if(init_rows!=init_cols){
		fprintf(stderr,"ERROR rows and cols must be equal\n");
		return return_error_val();
	}
	if(init_rows<=0){
		fprintf(stderr,"ERROR both rows and columns must have at least size 1\n");
		return return_error_val();
	}
	if(init_rows%num_proc!=0){
		fprintf(stderr,"ERROR num of processors must divide rows and columns up evenly.\n");
		fprintf(stderr,"Suggestions include rows and columns of size %d or %d\n",
				(init_rows/num_proc)*num_proc,(init_rows/num_proc)*(num_proc+1));
		return return_error_val();
	}
	if(flag_file==false && init_test==0){
		fprintf(stderr,"ERROR no input file has been specified if not running"
				" a test case must specify an input file.\n");
		return return_error_val();
	}
	if(POPT_Ret < -1){
		fprintf(stderr,"%s: %s\n",
				poptBadOption(POPT_Context,POPT_BADOPTION_NOALIAS),poptStrerror(POPT_Ret));
		return 1;
	} 
	poptFreeContext(POPT_Context);
	return 0;
}

/*------------------------------------------------------------------
 * Function:    print_help_message
 * Purpose:     Display help message to the screen
 */

void print_help_message(void){

    printf("\n");
		printf("To run the code use the following format\n");
		printf("\n");
		printf("mpiexec -np [num] ./executable [option] [arg]\n");
		printf("\n");
		printf("Where [num] is the number of processors\n");
		printf("\n");
		printf("Where [option] can be any one of the following or all of them:\n");
		printf("\n");
		printf(" -h     --help         To display help message no\n"
				   "                       argument is required\n");
		printf(" -v     --verbose      To display additional messages\n"
				   "                       during run no argument is required\n");
		printf(" -r     --rows         Specify the number of rows\n");
		printf(" -c     --cols         Specify the number of cols\n");
		printf(" -b     --block        Use the block implementation,\n"
				   "                       world divided up into panels\n");
		printf(" -g     --checkerboard Use checkerboard implementation,\n"
				   "                       world divided up into checkerboard\n"
				   "                       pattern\n");
		printf(" -s     --synchronous  Use synchronous communication\n");
		printf(" -a     --asynchronous Use asynchronous communication\n");
		printf(" -i     --iterations   Specify the number of iterations\n");
		printf(" -p     --print        Specifiy how often to print the \n"
				"                       bug sum\n");
		printf(" -d     --data         Print data file for matlab processing\n");
		printf(" -m     --print-pgm    Specify how often to print a pgm \n"
				"                       file.\n");
		printf(" -e     --ghost_count  Number of ghost rows and cols to pad \n"
				"                       the core matrix with.\n");
		printf(" -t     --test         Run a test scenario \n");
		printf(" -f     --file         Specify the input file name must \n"
				"                       a .pgm file.\n");

		printf("\n");
		printf("These options are optional default values exist for each\n"
				"parameter witht he exception of the file name:\n");
		printf(" * verbose      = false\n");
		printf(" * block        = true\n");
		printf(" * checkerboard = false\n");
		printf(" * synchronouse = false\n");
		printf(" * asynchronous = true\n");
		printf(" * iterations   = 10\n");
		printf(" * print        = 0\n");
		printf(" * data         = Not set\n");
		printf(" * print-pgm    = Not set\n");
		printf(" * ghost_count  = 1\n");
		printf(" * test         = 0\n");
		printf(" * file         = Not set\n");

		printf("\nNote argument format the brackets are not meant to be\n"
				"included in when running the executable but the quotation\n"
				"marks where indicated are necessary\n");
		printf(" --interations  [interger]\n");
		printf(" --print        [interger]\n");
		printf(" --print-pgm    [\"interger1-interger2\"]\n");
		printf(" --print-pgm    [\"interger1,interger2,interger3,etc\"]\n");
		printf(" --print-pgm    [\"interger\"]\n");
		printf(" --file         [filename]\n");

		printf("\n");
		printf("\nThis code is designed to execute Conway's Game of Life in parallel or serial.\n");
		printf("The following values can be set to choose the size of the\n"
				"matrix \n");
		printf(" * rows       - Total number of rows\n");
		printf(" * cols       - Total number of cols\n");

		printf("\nIf the world is being loaded from a .pgm file it is not\n"
           "necessary to specify the rows or cols because it is  \n"
           "contained in the .pgm file\n");


}

bool read_pgm(matrix * mat_a, matrix * mat_b){

//	pp_set_banner("pgm:readpgm");

	/* Open the file */
	if(my_rank==0){
		printf("Opening file %s\n",init_file);
		fflush(NULL);
	}
	FILE *fp = fopen( init_file, "r");
	
	if( !fp ){
		fprintf(stderr,"ERROR: The file '%s' cound not be opened.\n",init_file);
		return false;
	}

	if(*mat_a!=NULL){
		fprintf(stderr,"ERROR mat_a should be initialized to NULL before"
                       " passing to readpgm\n");
		return false;
	}
	
	if(*mat_b!=NULL){
		fprintf(stderr,"ERROR mat_b should be initialized to NULL before"
                       " passing to readpgm\n");
		return false;
	}
	/* Read the PGM header, which looks like this:
	 *  |P5        magic version number
	 *  |900 900       width height
	 *  |255         depth
	 */	
	char header[10];
	int width, height, depth;
	int rv = fscanf(fp, "%6s\n%i %i\n%i\n",header, &width, &height, &depth);
	if( rv!= 4){
		if(my_rank==0){
			fprintf(stderr,"ERROR: The file '%s' did not have a valid"
                     " PGM header\n",init_file);
		}
		return false;
	}
	if( my_rank==0 ){
		printf("%s: %s %i %i %i \n",init_file,header,width,height,depth);
		fflush(NULL);
	}
	/* Make sure the header is valid */
	if( strcmp( header, "P5") ){
		if(my_rank==0){
			fprintf(stderr,"ERROR: PGM file is not a valid P5 pixmap.\n");
		}
		return false;
	}

	if( depth != 255 ){
		if(my_rank==0){
			fprintf(stderr,"ERROR: PGM file has depth=%i, require depth=255.\n",
      						   depth);
		}
		return false;  
	}

	/* Make sure that the width and height are divisible by the number of 
     * processors in x and y directions */
	printf("width %d grid_c %d width mod grid_c %d\n",width,grid_dim_col,width%grid_dim_col);	

	if( width%grid_dim_col ){

		fprintf(stderr, "ERROR: %i pixel width cannot be divied into %i"
				" cols\n", width, grid_dim_col );
		return_error_val();
		return false;
	}
	printf("height %d grid_r %d height mod grid_r %d\n",height,grid_dim_row,height%grid_dim_row);
	fflush(NULL);	
	if( height%grid_dim_row ){
		fprintf(stderr, "ERROR: %d pixel height cannot  be divided into"
				" %d rows\n",	height, grid_dim_row );
		return_error_val();
		return false;
	}

	/* Divide the total image among the local processors */
	sub_matrix_dim_columns  = width  / grid_dim_col;
	sub_matrix_dim_rows     = height / grid_dim_row;

	/* Find out where my starting range is */
	int start_x = sub_matrix_dim_columns * my_grid_col;
	int start_y = sub_matrix_dim_rows    * my_grid_row;
	printf("Rank %d sub col %d my col %d sub row %d my row %d\n",
            my_rank,sub_matrix_dim_columns, my_grid_col, 
            sub_matrix_dim_rows, my_grid_row);	
	printf("Rank %d start x %d start y %d\n",my_rank,start_x, start_y);
	fflush(NULL);

	*mat_a = define_matrix_blockchecker();	
	*mat_b = define_matrix_blockchecker();	
	/* Read the data from the file. Save the local data to the local matrix. */
	int b, ll, lx, ly, y, x;
	for( y=0; y<height; y++){
		for( x=0; x<width; x++){
			/* Read the next character */
			b = fgetc(fp);
			if( b==EOF ){
				//pprintf("ERROR: Encountered EOF at [%i,%i]\n",y,x);
				return false;
			}
				
			/* From the PGM, black cells (b=0) are bugs, all other 
             * cells are background */
			if(b==0){
				b=1;
			}else{
				b=0;	
			}

			/* If the character is local then save it */

			if( x>=start_x && x< start_x + sub_matrix_dim_columns &&
                y>=start_y && y< start_y + sub_matrix_dim_rows    ){
			if(x==(width-1)){
//				printf("Rank %d x %d y %d start_x %d start_y %d sub col %d sub row %d\n",
//					 my_rank,x, y, start_x, start_y, sub_matrix_dim_columns, sub_matrix_dim_rows);
			}
				/* Calculate the local pixels */
				lx = x - start_x ;
				ly = y - start_y ;
				//printf("Rank %d lx %d ly %d\n",my_rank,lx,ly);	
				set_elem_core_matrix(*mat_a,ly,lx,b);
				set_elem_core_matrix(*mat_b,ly,lx,b);
			} /* Save local point */
			
		} /* For x */
	} /* for y */

	fclose(fp);
	//pp_reset_banner();
	return true;
}


/* Global constants that need to be initialized
 * before calling this function include 
 *  o my_grid_row
 *  o my_grid_col
 *  o grid_dim_row
 *  o grid_dim_col
 *  o my_rank
 *  o num_proc
 */
bool write_pgm(const_matrix mat_a, int count){


	char file_out[50];
	char header[50];
	int depth = 255;

	/* Write the data to the file */
	int b, ll, lx, ly, y, x;
	/* Create file name */
	printf("Rank %d\n",my_rank);
	sprintf(file_out,"fileout%d.pgm",count);

	/* Define width and height and finally the file header */
	int width = sub_matrix_dim_columns*grid_dim_col; 
	int height = sub_matrix_dim_rows*grid_dim_row;
	sprintf(header,"P5\n%d %d\n%d\n",width,height,depth);

	if(mat_a==NULL){
		fprintf(stderr,"ERROR mat_a should not be NULL when using write_pgm\n");
		return false;
	}
	
	#ifdef _MPI_H_
	
	/* Transfer core matrix data to char_array for binary format */
	int rows = get_rows_core_matrix(mat_a);
	int cols = get_cols_core_matrix(mat_a);
	/* Create temp char array */
	char char_array[rows*cols];
	for(x=0;x<rows;x++){
		for(y=0;y<cols;y++){
			if(get_elem_core_matrix(mat_a,x,y)==0.0){
				char_array[x*cols+y] = 255;
			}else{	
				char_array[x*cols+y] = 0;
			}
		}
	}
	/* Open the file */
	MPI_File fh;
	MPI_File_open(MPI_COMM_WORLD, file_out, 
                  MPI_MODE_CREATE | MPI_MODE_WRONLY,
                  MPI_INFO_NULL, &fh);
	/* Write the header to the file */
	if(my_rank==0){
		MPI_File_write(fh,(void*)header,strlen(header),MPI_CHAR,&status);
	}
	MPI_File_close(&fh);
	
	/* Define the MPI file type and etype */
	/* Define file type for the block implementation 
     * Makes sense in this case to write the file in blocks
     * that are equal to the row sizes */
	MPI_Datatype filetype;
	
	printf("Entering block on\n");

	if(block_on){

		/* With the darray approach */
		int array_of_gsizes[2];
		int array_of_distribs[2];
		int array_of_dargs[2];
		int array_of_psizes[2];
		int ndims = 2;

		/* Define the total dimensions */
		array_of_gsizes[0] = width*height;
		array_of_gsizes[1] = 1;

		/* Use the block decomposition */
		array_of_distribs[0] = MPI_DISTRIBUTE_BLOCK;
		array_of_distribs[1] = MPI_DISTRIBUTE_BLOCK;

		/* Set the number of dimensions each process will use */
		array_of_dargs[0] = rows*cols;
		array_of_dargs[1] = MPI_DISTRIBUTE_DFLT_DARG;

		/* Number of processors in each dimension */
		array_of_psizes[0] = num_proc;
		array_of_psizes[1] = 1;

		printf("gsizes[0] %d gsizes[1] %d\n"
               "dargs[0]  %d \n"
               "psizes[0] %d \n"
			   "strlen(header) %d \n",
               array_of_gsizes[0],
               array_of_gsizes[1],
               array_of_dargs[0],
               array_of_psizes[0],
			   (int)strlen(header));

		MPI_Type_create_darray(num_proc,
				my_rank,
				ndims,
				array_of_gsizes,
				array_of_distribs,
				array_of_dargs,
				array_of_psizes,
				MPI_ORDER_C,
				MPI_CHAR,
				&filetype);	

		MPI_Type_commit(&filetype);

		MPI_Offset displacement = strlen(header);
		MPI_File_open(MPI_COMM_WORLD, file_out, 
				MPI_MODE_CREATE | MPI_MODE_WRONLY,
				MPI_INFO_NULL, &fh);

		MPI_File_set_view(fh,displacement,MPI_CHAR,filetype,"native",MPI_INFO_NULL);
		
		MPI_File_write_all(fh,char_array,rows*cols,MPI_CHAR,&status);
	
		MPI_File_close(&fh);
	}else{
		/* With Checkerboard option */

		/* With the darray approach */
		int array_of_gsizes[2];
		int array_of_distribs[2];
		int array_of_dargs[2];
		int array_of_psizes[2];
		int ndims = 2;

		/* Define the total dimensions */
		array_of_gsizes[0] = width;
		array_of_gsizes[1] = height;

		/* Use the block decomposition */
		array_of_distribs[0] = MPI_DISTRIBUTE_BLOCK;
		array_of_distribs[1] = MPI_DISTRIBUTE_BLOCK;

		/* Set the number of dimensions each process will use */
		array_of_dargs[0] = rows;
		array_of_dargs[1] = cols;

		/* Number of processors in each dimension */
		array_of_psizes[0] = grid_dim_row;
		array_of_psizes[1] = grid_dim_col;

		printf("gsizes[0] %d gsizes[1] %d\n"
               "dargs[0]  %d \n"
               "psizes[0] %d \n"
			   "strlen(header) %d \n",
               array_of_gsizes[0],
               array_of_gsizes[1],
               array_of_dargs[0],
               array_of_psizes[0],
			   (int)strlen(header));

		MPI_Type_create_darray(num_proc,
				my_rank,
				ndims,
				array_of_gsizes,
				array_of_distribs,
				array_of_dargs,
				array_of_psizes,
				MPI_ORDER_C,
				MPI_CHAR,
				&filetype);	

		MPI_Type_commit(&filetype);

		MPI_Offset displacement = strlen(header);
		MPI_File_open(MPI_COMM_WORLD, file_out, 
				MPI_MODE_CREATE | MPI_MODE_WRONLY,
				MPI_INFO_NULL, &fh);

		MPI_File_set_view(fh,displacement,MPI_CHAR,filetype,
                          "native",MPI_INFO_NULL);
		
		MPI_File_write_all(fh,char_array,rows*cols,MPI_CHAR,&status);
	
		MPI_File_close(&fh);

	}
	#else
	
	FILE *fp = fopen( file_out, "w");
	if( !fp ){
		//pprintf( "Error: The file '%s' could not be opened.\n", file_out );
		return false;
	}
	fprintf(fp,"%s",header);
	/* Find out where my starting range is */
	int start_x = sub_matrix_dim_columns * my_grid_col;
	int start_y = sub_matrix_dim_rows    * my_grid_row;

	for( y=0; y<height; y++){
		for( x=0; x<width; x++){
			/* Read the next character */
	
			/* If the character is local then save it */
			if( x>=start_x && x< start_x + sub_matrix_dim_columns &&
                y>=start_y && y< start_y + sub_matrix_dim_rows    ){
			
				/* Calculate the local pixels */
				lx = x - start_x ;
				ly = y - start_y ;
			
				/* From the PGM, black cells (b=0) are bugs, all other 
				 * cells are background */
				if(get_elem_core_matrix(mat_a,ly,lx)==1){
					fputc(0,fp);
				}else if(get_elem_core_matrix(mat_a,ly,lx)==0){
					fputc(255,fp);
				}else{
					printf("ERROR Matrix is only meant to contain 0s and 1s\n");
					return false;
				}	
			} /* Save local point */
		} /* For x */
	} /* for y */
	fclose(fp);
	#endif

	return true;
}

matrix define_matrix_blockchecker(void){

	int above_ghost;
	int below_ghost;

	/* For the Block Panel decompostion */
	if(block_on){
		if( my_rank==0 && my_rank==(num_proc-1)){
			/* Single Processor */
			above_ghost = 0;
			below_ghost = 0;

		}else if(my_rank==0){
			/* rank 0 is at the top of the world thus we do not need a ghost 
			 * row for it's top row */
			above_ghost = 0;
			below_ghost = ghost_count;

		}else if(my_rank==(num_proc-1)){
			/* rank == (num_proc-1) is at the bottom of the world thus we can
			 * omit the bottom row */
			above_ghost = ghost_count;
			below_ghost = 0;

		}else{
			/* Panels exist above and below */
			above_ghost = ghost_count;
			below_ghost = ghost_count;
		}

		printf("sub_row %d sub_col %d ",sub_matrix_dim_rows,sub_matrix_dim_columns);
		printf("above_g %d below_g %d ",above_ghost,below_ghost);
		return new_matrix_row_panel( sub_matrix_dim_rows, 
				                     sub_matrix_dim_columns, 
				                     above_ghost, below_ghost);
	}else{
		/* For the checkerboard decomposition */
		int right_ghost;
		int left_ghost;	

		/* Predefine the ghost_rows to exist */
		left_ghost  = ghost_count;
		right_ghost = ghost_count;
		above_ghost = ghost_count;
		below_ghost = ghost_count;
	
		/* Determine if my_rank is on the edges  */
		if(my_grid_col==0){
			left_ghost = 0;		
		}	
		if(my_grid_col==(grid_dim_col-1)){
			right_ghost = 0;
		}
		if(my_grid_row==0){
			above_ghost = 0;
		}
		if(my_grid_row==(grid_dim_row-1)){
			below_ghost = 0;
		}
	
		return new_matrix_checkerboard_elem(sub_matrix_dim_rows,
                                            sub_matrix_dim_columns,
                                            above_ghost, below_ghost,
                                            right_ghost, left_ghost );
	}

}

int calculate_statistics(int iteration){

	if(num_proc>2){
		MPI_Reduce(&elapsed,&max_t, 1, MPI_DOUBLE,
                   MPI_MAX, 0, MPI_COMM_WORLD);
	}else{
		max_t = elapsed;
	}
	if(my_rank==0){
		stat_sum += max_t;
		running_mean = stat_sum/(iteration+1);
		standard_dev = standard_deviation(running_mean, iteration);
		margin_err = t_val*standard_dev/pow((double)iteration,1.0/2.0);	
	}
	return 0;
}

double standard_deviation(double mean, int iteration){
	return pow(pow(max_t-mean,2)/((double)iteration+1.0),(1.0/2.0));
}

int generate_performance_file( char * file_name, int iter, int bytes){
	
	FILE *fp;
	fp = fopen(file_name,"a+");
	fprintf(fp,"Iter %8d ",iter);
	fprintf(fp,"Bytes %8d ",bytes);
	fprintf(fp,"Mean %12g ",running_mean);
	fprintf(fp,"Dev %12g ",standard_dev);
	fprintf(fp,"Error+- %12g\n",margin_err);
	fclose(fp);
	return 0;
}
