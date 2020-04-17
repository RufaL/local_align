 /*
  * Generate random strings of length 'L' with just the characters A, G, C and T
  * For L = 100, 150 and 300
  */
 #include "stdio.h"
 #include "string.h"
 #include "stdlib.h"
 #include "stdint.h"
 #include "time.h"

 #define L 100
 #define n_seq 10

 int main(){
 int lower = 1, upper = 4, count, num;
 char seq1[300], seq2[300];
 FILE *output1, *output2;
 output1 = fopen("seq1_out.txt","wb");
 output2 = fopen("seq2_out.txt","wb");
 count  = 2*L;
 for(int i=0; i<n_seq; i++){
         srand(i);
	 for(int j=0; j<count; j++){
	   num = (rand() % (upper-lower+1)) + lower;
           //printf("Random no.:%d\n",num);
	   switch(num){
	   	case 1: if(j < L)
			  seq1[j] = 'A';
                        else
                          seq2[j-L] = 'A';
	   	        break;
	   	case 2: if(j < L)
			  seq1[j] = 'T';
                        else
                          seq2[j-L] = 'T';
	   	        break;
	   	case 3: if(j < L)
			  seq1[j] = 'G';
                        else
                          seq2[j-L] = 'G';
	   	        break;
	   	case 4: if(j < L)
			  seq1[j] = 'C';
                        else
                          seq2[j-L] = 'C';
	   	        break;
	   }
	 }
	fwrite(seq1, sizeof(char), L, output1);
        fprintf(output1,"\n");
        fwrite(seq2, sizeof(char), L, output2);
        fprintf(output2,"\n");
	
 }	 
 fclose(output1);
 fclose(output2);



 return 0;
 }
