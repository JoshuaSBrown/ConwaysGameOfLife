#include <stdio.h>

int main(void){

	FILE * fp;
	fp = fopen("test3.pgm","w");
	unsigned char ch = 0xc5;

	fprintf(fp,"P5\n%d %d\n%d\n",900,900,255);
	fputc(255,fp);
	fputc(255,fp);
	fputc(255,fp);
	fputc(255,fp);
	fputc(255,fp);
	fputc(255,fp);
	fputc(255,fp);
	fputc(255,fp);
	fputc(255,fp);
	fputc(255,fp);
	fputc(255,fp);
	fputc(255,fp);
	fputc(255,fp);
	fputc(255,fp);
	fputc(255,fp);
	fputc(255,fp);

	int i = 1;
	fputc(i,fp);
/*	fprintf(fp,"%u\n",'a'-97);
	fprintf(fp,"%u\n",'a'-97);
	fprintf(fp,"%u\n",'a'-97);
	fprintf(fp,"%u\n",'a'-97);
	fprintf(fp,"%u\n",'a'+158);
	fprintf(fp,"%u\n",'a'+158);
	fprintf(fp,"%u\n",'a'+158);
	fprintf(fp,"%u\n",'a'+158);
*/
	printf("\n");
	printf("%u", 'a'+158);
	printf("%u", 'a'+158);
	printf("%u", 'a'+158);

	printf("%u", 'a'-97);
	printf("%u", 'a'-97);
	printf("%u", 'a'-97);
	printf("\n");
	fclose(fp);

	return 0;
}
