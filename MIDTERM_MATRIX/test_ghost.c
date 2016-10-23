#include <stdio.h>
#include <stdlib.h>
#include "midterm_matrix.h"

float add_all(matrix * orig_mat, int i,int j){
	
	float sum = 0.0;
	float value;
	float NW;
	float N;
	float NE;
	float E;
	float W;
	float SW;
	float S;
	float SE;

	value = get_elem(*orig_mat,i,j);
	if(within_matrix(*orig_mat,i-1,j-1)){
		NW = get_elem(*orig_mat,i-1,j-1);
		sum += NW;
	}
	if(within_matrix(*orig_mat,i,j-1)){
		N = get_elem(*orig_mat,i,j-1);
		sum += N;
	}
	if(within_matrix(*orig_mat,i+1,j-1)){
		NE = get_elem(*orig_mat,i+1,j-1);
		sum += NE;
	}
	if(within_matrix(*orig_mat,i-1,j)){
		W = get_elem(*orig_mat,i-1,j);
		sum += W;
	}
	if(within_matrix(*orig_mat,i+1,j)){
		E = get_elem(*orig_mat,i+1,j);
		sum += E;
	}
	if(within_matrix(*orig_mat,i-1,j+1)){
		SW = get_elem(*orig_mat,i-1,j+1);
		sum += SW;
	}
	if(within_matrix(*orig_mat,i,j+1)){
		S = get_elem(*orig_mat,i,j+1);
		sum += S;
	}
	if(within_matrix(*orig_mat,i+1,j+1)){
		SE = get_elem(*orig_mat,i+1,j+1);
		sum += SE;
	}

	return sum;
}


int main(void){
	
	int ghost = 2;

	matrix orig = new_matrix_checkerboard_elem(12,12,0,ghost,0,ghost);
	matrix mat1 = new_matrix_checkerboard_elem(12,12,0,ghost,0,ghost);
	matrix mat2 = new_matrix_checkerboard_elem(12,12,0,ghost,0,ghost);
	set_all_elems(orig,1.0);
	set_all_elems(mat1,1.0);
	set_all_elems(mat2,1.0);
	print_matrix_compressed(mat1);
	matrix temp1 = mat1;
	matrix temp2 = mat2;
	matrix new;
	int ghost_left  = get_cols_left_ghost_matrix(mat1);
	int ghost_right = get_cols_right_ghost_matrix(mat1);
	int ghost_above = get_rows_above_ghost_matrix(mat1);
	int ghost_below = get_rows_below_ghost_matrix(mat1);

	int col_start;
	int col_end;
	int row_start;
	int row_end;

	float val;

	int i, j, k;
	for(k=0;k<ghost;k++){


		if(ghost_below){
			row_end   = get_row_end(temp1)-(ghost_below+k);
		}else{
			row_end = get_rows(temp1)-1;
		}
		if(ghost_above){
			row_start = get_row_start(temp1)+(ghost_above-1+k);
		}else{
			row_start = 0;
		}
		col_start = get_col_start(temp1)-(ghost_left-k-1);
		col_end   = get_col_start(temp1)+(ghost_left*2-k-1);
		printf("cs %d ce %d rs %d re %d\n",
				col_start,col_end,row_start,row_end);
		for(j=col_start;j<col_end;j++){
			for(i=row_start;i<row_end+1;i++){
				//set_elem(temp1,i,j,(float)k+1.0);		
				//val = get_elem(temp1,i,j);
				val = add_all(&temp1,i,j);
				set_elem(temp2,i,j,val);
			}
		}		
		print_matrix_compressed(temp2);
		getchar();

		if(ghost_below){
			row_end   = get_row_end(temp1)-(ghost_below+k);
		}else{
			row_end = get_rows(temp1)-1;
		}
		if(ghost_above){
			row_start = get_row_start(temp1)+(ghost_above-1+k);
		}else{
			row_start = 0;
		}	
		col_start = get_col_end(temp1)-(ghost_right*2-k-1);
		col_end   = get_col_end(temp1)+(ghost_right-k-1);
		for(j=col_start;j<col_end;j++){
			for(i=row_start;i<row_end+1;i++){
				//set_elem(temp1,i,j,(float)k+4.0);
				//val = get_elem(temp1,i,j);
				val = add_all(&temp1,i,j);
				set_elem(temp2,i,j,val);
			}
		}

		print_matrix_compressed(temp1);
		print_matrix_compressed(temp2);
		getchar();
		/* Top Row */
		if(ghost_above){
			if(ghost_left==0){
				col_start=0;
			}else{
				col_start = get_col_start(temp1)-(ghost_left-k-1);
			}
			if(ghost_right==0){
				col_end = get_cols(temp1);
			}else{
				col_end   = get_col_end(temp1)+(ghost_right-k-1);
			}
			row_start = get_row_start(temp1)-(ghost_above-k-1);
			row_end   = get_row_start(temp1)+(ghost_above-k-2);
			for(j=col_start;j<col_end;j++){
				for(i=row_start;i<row_end+1;i++){
					//set_elem(temp1,i,j,(float)k+7.0);
					//val = get_elem(temp1,i,j);
					val = add_all(&temp1,i,j);
					set_elem(temp2,i,j,val);
				}
			}
			printf("one\n");
			if(ghost_left==0){
				col_start=0;
			}else{
				col_start = get_col_start(temp1)+(ghost_left*2-k-1);
			}
			if(ghost_right==0){
				col_end = get_cols(temp1);
			}else{
				col_end   = get_col_end(temp1)-(ghost_right*2-k-1);
			}
			print_matrix_compressed(temp2);
			print_matrix_compressed(temp1);
			printf("col_start %d col_end %d\n",col_start,col_end);
			getchar();
			row_start = row_end+1;
			row_end   = row_end+ghost_above;
			for(j=col_start;j<col_end;j++){
				for(i=row_start;i<row_end+1;i++){
					//set_elem(temp1,i,j,(float)k+16.0);
					//val = get_elem(temp1,i,j);
					val = add_all(&temp1,i,j);
					set_elem(temp2,i,j,val);
				}
			}
			printf("two\n");
			print_matrix_compressed(temp2);
			getchar();
		}
		/* Bottom Row */
		if(ghost_below){
			if(ghost_left==0){
				col_start=0;
			}else{
				col_start = get_col_start(temp1)-(ghost_left-k-1);
			}
			if(ghost_right==0){
				col_end = get_cols(temp1);
			}else{
				col_end   = get_col_end(temp1)+(ghost_right-k-1);
			}

			row_start = get_row_end(temp1)-(ghost_below-k-1);
			row_end   = get_row_end(temp1)+(ghost_below-k-2);
			printf("cs %d ce %d rs %d re %d\n",
					col_start,col_end,row_start,row_end);
			for(j=col_start;j<col_end;j++){
				for(i=row_start;i<row_end+1;i++){
					//set_elem(temp1,i,j,(float)k+10.0);
					//val = get_elem(temp1,i,j);
					val = add_all(&temp1,i,j);
					set_elem(temp2,i,j,val);
				}
			}
			print_matrix_compressed(temp2);
			getchar();
			
			if(ghost_left==0){
				col_start=0;
			}else{
				col_start = get_col_start(temp1)+(ghost_left*2-k-1);
			}
			if(ghost_right==0){
				col_end = get_cols(temp1);
			}else{
				col_end   = get_col_end(temp1)-(ghost_right*2-k-1);
			}
			row_end   = row_start;
			row_start = row_start-ghost_below;
			for(j=col_start;j<col_end;j++){
				for(i=row_start;i<row_end;i++){
					//set_elem(temp1,i,j,(float)k+13.0);
					//val = get_elem(temp1,i,j);
					val = add_all(&temp1,i,j);
					set_elem(temp2,i,j,val);
				}
			}

			print_matrix_compressed(temp2);
			getchar();
		}

		/* Top Corners */
		if(ghost_above){
			/* North West Corner */		
			row_start = get_row_start(temp1)-k+ghost_above-1;
			row_end   = get_row_start(temp1)+k+ghost_above-1;
			col_start = get_col_start(temp1)-(ghost_left-k-1);
			col_end   = get_col_start(temp1)+(ghost_left*2-k-1);

			for(j=col_start;j<col_end;j++){
				for(i=row_start;i<row_end;i++){
					//set_elem(temp1,i,j,(float)k+20.0);			
					//val = get_elem(temp1,i,j);
					val = add_all(&temp1,i,j);
					set_elem(temp2,i,j,val);
				}
			}		
			print_matrix_compressed(temp2);
			getchar();
			
			/* North East Corner */
			row_start = get_row_start(temp1)-k+ghost_above-1;
			row_end   = get_row_start(temp1)+k+ghost_above-1;
			col_start = get_col_end(temp1)-(ghost_right*2-k-1);
			col_end   = get_col_end(temp1)+(ghost_right-k-1);

			for(j=col_start;j<col_end;j++){
				for(i=row_start;i<row_end;i++){
					//set_elem(temp1,i,j,(float)k+20.0);			
					//val = get_elem(temp1,i,j);
					val = add_all(&temp1,i,j);
					set_elem(temp2,i,j,val);
				}
			}		
			print_matrix_compressed(temp2);
			getchar();
		}

		if(ghost_below){
			/* Sourth East Corner */
			row_start = get_row_end(temp1)-k-(ghost_below-1);
			row_end   = get_row_end(temp1)+k-(ghost_below-1);
			col_start = get_col_end(temp1)-(ghost_right*2-k-1);
			col_end   = get_col_end(temp1)+(ghost_right-k-1);

			for(j=col_start;j<col_end;j++){
				for(i=row_start;i<row_end;i++){
					//set_elem(temp1,i,j,(float)k+25.0);			
					//val = get_elem(temp1,i,j);
					val = add_all(&temp1,i,j);
					set_elem(temp2,i,j,val);
				}
			}		
			print_matrix_compressed(temp2);
			getchar();

			/* South West Corner */
			row_start = get_row_end(temp1)-k-(ghost_below-1);
			row_end   = get_row_end(temp1)+k-(ghost_below-1);
			col_start = get_col_start(temp1)-(ghost_left-k-1);
			col_end   = get_col_start(temp1)+(ghost_left*2-k-1);

			for(j=col_start;j<col_end;j++){
				for(i=row_start;i<row_end;i++){
					//set_elem(temp1,i,j,(float)k+25.0);		
					//val = get_elem(temp1,i,j);
					val = add_all(&temp1,i,j);
					set_elem(temp2,i,j,val);

				}
			}		
			print_matrix_compressed(temp2);
			getchar();
		}

		if(k%2==0){
			temp1 = mat2;
			temp2 = mat1;
		}else{
			temp1 = mat1;
			temp2 = mat2;
		}
	}

	/* At this point all the innter columns and rows 
     * have been updated and the ghost rows and columns
     * can be updated*/
	printf("Calculating Inner\n");
	temp1 = orig;
	if(ghost==1){
		new = mat2;
		temp2 = mat1;
	}else{
		if(ghost%2==0){
			temp2 = mat2;
			new = mat1;
		}else{
			temp2 = mat1;
			new = mat2;
		}
	}


	if(ghost%2==1){
		/* Odd case */
		for(k=0;k<(ghost/2+1);k++){
			if(ghost_above==0){
				row_start = 0;
			}else{
				row_start = get_row_start(orig)+1+k;
			}
			if(ghost_below==0){
				row_end = get_rows(orig);
			}else{
				row_end = get_row_end(orig)-1-k;
			}
			if(ghost_left==0){
				col_start=0;
			}else{
				col_start = get_col_start(orig)+1+k;
			}
			if(ghost_right==0){
				col_end = get_cols(orig);
			}else{
				col_end = get_col_end(orig)-1-k;
			}

			if(ghost/2==k){
				for(j=col_start;j<col_end/2;j++){
					for(i=row_start;i<row_end;i++){
						//val = get_elem(temp1,i,j);
						val = add_all(&temp1,i,j);
						if(k==0 && (ghost)!=1){
							set_elem(temp2,i,j,val);
						}else if(k==(ghost-1)){
							set_elem(new,i,j,val);	
						}else{
							set_elem(temp2,i,j,val);	
						}
					}
				}
			}else{
				for(j=col_start;j<col_end;j++){

					for(i=row_start;i<row_end;i++){
						//val = get_elem(temp1,i,j);
						val = add_all(&temp1,i,j);
						if(k==0 && (ghost)!=1){
							set_elem(temp2,i,j,val);
						}else if(k==(ghost-1)){
							set_elem(new,i,j,val);	
						}else{
							set_elem(temp2,i,j,val);	
						}
					}
				}
				if(k%2==0){
					temp1 = temp2;
					temp2 = orig;
				}else{
					temp2 = temp1;
					temp1 = orig;
				}
			}
		}

		/* Halfway point */
		for(k=ghost/2;k<ghost;k++){
		/* Odd */
			if(ghost_above==0){
				row_start = 0;
			}else{
				row_start = get_row_start(orig)+1+k;
			}
			if(ghost_below==0){
				row_end = get_rows(orig);
			}else{
				row_end = get_row_end(orig)-1-k;
			}
			if(ghost_left==0){
				col_start=0;
			}else{
				col_start = get_col_start(orig)+1+k;
			}
			if(ghost_right==0){
				col_end = get_cols(orig);
			}else{
				col_end = get_col_end(orig)-1-k;
			}

			if(ghost/2==k){
				for(j=col_end/2;j<col_end;j++){
					for(i=row_start;i<row_end;i++){
						//val = get_elem(temp1,i,j);
						val = add_all(&temp1,i,j);
						if(k==0 && (ghost)!=1){
							set_elem(temp2,i,j,val);
						}else if(k==(ghost-1)){
							set_elem(new,i,j,val);	
						}else{
							set_elem(temp2,i,j,val);	
						}
					}
				}
			}else{
				for(j=col_start;j<col_end;j++){

					for(i=row_start;i<row_end;i++){
						//val = get_elem(temp1,i,j);
						val = add_all(&temp1,i,j);
						if(k==0 && (ghost)!=1){
							set_elem(temp2,i,j,val);
						}else if(k==(ghost-1)){
							set_elem(new,i,j,val);	
						}else{
							set_elem(temp2,i,j,val);	
						}
					}
				}

			}
			print_matrix_compressed(temp2);
			getchar();
			if(k%2==0){
				temp1 = temp2;
				temp2 = orig;
			}else{
				temp2 = temp1;
				temp1 = orig;
			}

		}
	}else{
		/* Even case */
		printf("Even ghost rows \n");

		for(k=0;k<(ghost/2+1);k++){
			if(ghost_above==0){
				row_start = 0;
			}else{
				row_start = get_row_start(orig)+1+k;
			}
			if(ghost_below==0){
				row_end = get_rows(orig);
			}else{
				row_end = get_row_end(orig)-1-k;
			}
			if(ghost_left==0){
				col_start=0;
			}else{
				col_start = get_col_start(orig)+1+k;
			}
			if(ghost_right==0){
				col_end = get_cols(orig);
			}else{
				col_end = get_col_end(orig)-1-k;
			}

			for(j=col_start;j<col_end;j++){

				for(i=row_start;i<row_end;i++){
					//val = get_elem(temp1,i,j);
					val = add_all(&temp1,i,j);
					if(k==0 && (ghost)!=1){
						set_elem(temp2,i,j,val);
					}else if(k==(ghost-1)){
						set_elem(new,i,j,val);	
					}else{
						set_elem(temp2,i,j,val);	
					}
				}
			}
		
			if(k%2==0){
				temp1 = temp2;
				temp2 = orig;
			}else{
				temp2 = temp1;
				temp1 = orig;
			}
		}

		printf("temp1\n");
		getchar();
		print_matrix_compressed(temp1);
		printf("orig\n");
		getchar();
		print_matrix_compressed(orig);
		/* Halfway point */
		for(k=(ghost/2+1);k<ghost;k++){
			if(ghost_above==0){
				row_start = 0;
			}else{
				row_start = get_row_start(orig)+1+k;
			}
			if(ghost_below==0){
				row_end = get_rows(orig);
			}else{
				row_end = get_row_end(orig)-1-k;
			}
			if(ghost_left==0){
				col_start=0;
			}else{
				col_start = get_col_start(orig)+1+k;
			}
			if(ghost_right==0){
				col_end = get_cols(orig);
			}else{
				col_end = get_col_end(orig)-1-k;
			}

			for(j=col_start;j<col_end;j++){

				for(i=row_start;i<row_end;i++){
					//val = get_elem(temp1,i,j);
					val = add_all(&temp1,i,j);
					if(k==0 && (ghost)!=1){
						set_elem(temp2,i,j,val);
					}else if(k==(ghost-1)){
						set_elem(new,i,j,val);	
					}else{
						set_elem(temp2,i,j,val);	
					}
				}
			}

			if(k%2==0){
				temp1 = temp2;
				temp2 = orig;
			}else{
				temp2 = temp1;
				temp1 = orig;
			}
		}
		print_matrix_compressed(temp2);
		getchar();

	}

	printf("ghost %d\n",ghost);

	print_matrix_compressed(new);
	delete_matrix(&mat1);
	delete_matrix(&mat2);
	delete_matrix(&orig);

	return 0;
}
