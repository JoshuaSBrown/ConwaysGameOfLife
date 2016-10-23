#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "../MIDTERM_MATRIX/midterm_matrix.h"
#include "../MIDTERM_CONSTANTS/midterm_constants.h"
#include "../MIDTERM_ERROR/midterm_error.h"

#include "mpi.h"


/* Will let the user know if is a single iteration we will print at
 * a list of values we will print out, or a range of iterations we
 * will look at will return:
 *   -1 Malformed string input 
 *    0 Single value input
 *    1 List comma deliminated
 *    2 range dash deliminated
 */
int single_list_range(char * str){

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
				return -1;
			}
		}
	}
	
	if(dash==1 && comma==1){
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
	return -1;
}

int string_check_for_dashs(char * str){

    if(!str){
        fprintf(stderr,"ERROR string is NULL in string_check_for_dashs\n");
        exit(1);
    }

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

    if(!str){
        fprintf(stderr,"ERROR str is NULL in string_check_for_commas\n");
        exit(1);
    }
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
        exit(1);
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

    if(!mat){
        fprintf(stderr,"ERROR mat is NULL in pgm_read_values\n");
        exit(1);
    }
    if(!str){
        fprintf(stderr,"ERROR str is NULL in pgm_read_values\n");
        exit(1);
    }

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

	        if(get_elem(mat,0,0)>get_elem(mat,1,0)){
            fprintf(stderr,"ERROR When entering a range the smaller number must come first\n");
            exit(1);
        }
        return 0;
    }

    /* Single interger no dashs or commas */
    set_elem(mat,0,0,(float)atoi(str));
	free(str_copy);
	free(str_copy2);
    return 0;
}
	


int main(void){

	int i;
	int rv;
	char * str1 = "Hello";
	
	rv = string_check_for_dashs(str1);
	assert(rv==0);
	rv = string_check_for_dashs("Hello,Oh Hi");
	assert(rv==0);
	rv = string_check_for_dashs("Hello-Oh Hi");
	assert(rv==1);
	rv = string_check_for_dashs("Hello-Oh-Hi");
	assert(rv==2);
	rv = string_check_for_commas("Hello Yo");
	assert(rv==0);
	rv = string_check_for_commas(",Hi,,How are you");
	assert(rv==3);
	
	int size = pgm_input_number_values("1423");
	assert(size==0);
	size = pgm_input_number_values("232,14");
	assert(size==1);
	size = pgm_input_number_values("4324, 12343, 5434");
	assert(size==2);
	size = pgm_input_number_values("435-13");
	assert(size==1);

	matrix mat1;
	char * str2 = "323,1345,595,23";
	size = pgm_input_number_values(str2);
	assert(size==3);
	mat1 = new_matrix(size+1,1);
	rv = pgm_read_values_to_array(mat1,str2);
	assert(rv==0);
	print_matrix(mat1);
	delete_matrix(&mat1);
	char * str3 = "434,13,134,";
	size = pgm_input_number_values(str3);
	assert(size==3);
	mat1 = new_matrix(size+1,1);
	rv = pgm_read_values_to_array(mat1,str3);
	print_matrix(mat1);
	delete_matrix(&mat1);	
	char * str4 = "134-500";
	size = pgm_input_number_values(str4);
	assert(size==1);
	mat1 = new_matrix(size+1,1);
	rv = pgm_read_values_to_array(mat1,str4);
	print_matrix(mat1);
	delete_matrix(&mat1);	
	char * str5 = "500";
	size = pgm_input_number_values(str5);
	assert(size==0);
	mat1 = new_matrix(size+1,1);
	rv = pgm_read_values_to_array(mat1,str5);
	print_matrix(mat1);
	delete_matrix(&mat1);	


	rv = single_list_range("");
	assert(rv==-1);
	rv = single_list_range("43");
	assert(rv==0);
	rv = single_list_range("43,13-23");
	assert(rv==-1);
	rv = single_list_range("235,1345");
	assert(rv==1);
	rv = single_list_range("134-1345");
	assert(rv==2);
	rv = single_list_range("2334-24566-234566");
	assert(rv==-1);
	//int array1[size+1];
	//rv = pgm_read_values_to_array(array1,str2);
	//assert(rv==0);
//	for(i=0;i<size+1;i++){
//		printf("%d",array1[i]);
//	}

	return 0;
}
