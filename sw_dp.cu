/*
 * The Smith-Waterman algorithm, is a dynamic programming algorithm were the DP matrices 
 * involved in the computation are calculated dynamically. There are 3 DP: M, X, and Y 
 * each contributing a score from one of the three directions an entry in the SW scoring 
 * matrix can obtain. With the SW algorithm we implement affine-gap penalty scoring, thus
 * working towards a local alignment algortihm with affine-gap penalty as in the seed 
 * extension stage of BWA-MEM sequencing algorithm. 
 */
 #include "stdio.h"
 #include "string.h"
 #include "stdlib.h"
 #include "stdint.h"
 #include "swalign.h"
 #include <iostream>
 #include "cuda_runtime.h"

/*const int for penalty*/
const int penalty = gap_open + gap_extn;    
#define size ((L+1)/8 +1)


__host__ __device__ void init_DP(int16_t M[][L+1], int16_t X[][L+1], int16_t Y[][L+1]){
	M[0][0] = 0;
	X[0][0] = -1000;
	Y[0][0] = -1000;
	for(int i=1; i <L+1; i++){
		M[i][0] = 0;
		X[i][0] = -1000;   //Just a large negative number
		Y[i][0] = -1000;
	}

	for(int j=1; j< L+1; j++){
		M[0][j] = 0;
		X[0][j] = -1000;   //Just a large negative number
		Y[0][j] = -1000;
	}
}

__host__ __device__ void unpacking(uint32_t *s1, char *seq1_out){
    uint8_t c;
    int m=0;
     
    for(int i=0; i<size; i++){
	    c = s1[i] & 0x0F;
        for(int j=0; (c!=0 && j<8); ){
          
           switch(c){
            case 0x1: seq1_out[m] = 'A';
                          break;
	    case 0x3: seq1_out[m] = 'C';
                         break;		      
            case 0x7: seq1_out[m] = 'G';
                          break;
            case 0x4: seq1_out[m] = 'T';
                          break;
            case 0xD: seq1_out[m] = '-';
                          break;
            case 0xE: seq1_out[m] = '.';
                          break;	      
            case 0xA: seq1_out[m] = '*';
                          break;
            }
	   m++;
	   j++;
	   c = (s1[i] >> 4*j) & 0x0F;
        
        }
    }
}

__global__ void read_align(char *sq1, char *sq2, char *seq1_out, char *seq2_out){
    
   int seq_i, sq_i;
   uint16_t Score_Matrix[L+1][L+1];
   uint8_t Dir[L][L/2];
   int16_t M[2][L+1], X[2][L+1], Y[2][L+1];  //DP matrices
   int A, B, S_I;
   uint32_t s1_out[size], s2_out[size];  
   uint32_t  seq1[size], seq2[size];

   int test;

   int index = blockIdx.x * blockDim.x +threadIdx.x;
   
   if(index < no_seq)
   {   
        seq_i = index * size;
	sq_i = index * (L+1);
	sq1[sq_i] = '-';
	sq2[sq_i] = '-';
        /*Start scoring*/
       
        init_DP(M, X, Y);
        /*data packing*/
	    int p, j=0;
	    seq1[0] = 0x0;
	    seq2[0] = 0x0;
	    for(int i=0; i<(L+1); i++){
		p = i%8;    
	        seq1[j] |= (sq1[i + sq_i] & 0x0F) << 4*p;		  
	        seq2[j] |= (sq2[i + sq_i] & 0x0F) << 4*p;
				 
		if(p==7){
		  ++j;
		  seq1[j] = 0;
		  seq2[j] = 0;
		} 
	    }
       
	 //unpacking(seq1, &seq1_out[sq_i]);
	 //unpacking(seq2, &seq2_out[sq_i]);
       

      /*Initializing score martix*/  
      
            Score_Matrix[0][0] = 0;
            for(int j=1; j<L+1; j++){
              Score_Matrix[0][j] = 0;
              Score_Matrix[j][0] = 0;
            }
       //A = M[0][0];
       //seq1_out[A] = 'Z';
	  
 /*Compute DP matrices */
      
    int16_t M_max =0, X_max, Y_max;
    int16_t M_x, M_y, M_m;
    int16_t match_score;
    int si, sj;
    int r, c;
    int sw_max;
    uint8_t e1, e2;
    Dir[0][0] = 0;
    int d_I=0, d_J=0, d_p=0;

    for(int I = 1; I < L+1; I=I+8){
       for(int J = 1; J <L+1; J=J+8){
	    r = I/8;
            c = J/8;   
            for(int i=0; i<8; i++){
		    if(i == 7)
			e1 = seq1[r+1] & 0x0F;
		    else 
			e1 = (seq1[r] >> 4*(i+1)) & 0x0F;
                for(int j=0; j<8; j++){
			if(j == 7)
			   e2 = seq2[c+1] & 0x0F;
			else 
			   e2 = (seq2[c] >> 4*(j+1)) & 0x0F;
		    if((I+i < L+1) && (J+j < L+1)){	
		         
                       if(e1 == e2)
                        match_score = match;
                       else
                        match_score = mismatch;
                           
                       M_m = M[(I+i-1)%2][J+j-1] + match_score;
                       M_x = X[(I+i-1)%2][J+j-1] + match_score;
                       M_y = Y[(I+i-1)%2][J+j-1] + match_score;

		        M_max = 0;
                        if(M_m >= M_x && M_m >= M_y && M_m > 0) 
                            M_max = M_m;
                        else if(M_x >= M_m && M_x >= M_y && M_x > 0)
                            M_max = M_x;
                             else if(M_y >= M_m && M_y >= M_x && M_y > 0)
                                 M_max = M_y;

                        M[(I+i)%2][J+j] =  M_max;
                         
                        Y_max = gap_extn + Y[(I+i)%2][J+j-1];
                        if(penalty + M[(I+i)%2][J+j-1] > Y_max)
                        Y_max = M[(I+i)%2][J+j-1] + penalty;

                        Y[(I+i)%2][J+j] = Y_max;

                        X_max = gap_extn + X[(I+i-1)%2][J+j];
                        if(penalty + M[(I+i-1)%2][J+j] > X_max)
                        X_max = M[(I+i-1)%2][J+j] + penalty;

                        X[(I+i)%2][J+j] = X_max;

                         
                        if(X_max >= Y_max && X_max >= M_max){
                          Score_Matrix[I+i][J+j] = X_max;
                          Dir[d_I][d_J] |= (0x02 << 4*d_p);
                        }
                        else if(Y_max >= X_max && Y_max >= M_max){
                            //Score_Matrix[I+i][J+j] = Y_max;
                            Dir[d_I][d_J] |= (0x03 << 4*d_p);
                             }
                         else if(M_max >= X_max && M_max >= Y_max){
                             //Score_Matrix[I+i][J+j] = M_max;
                             Dir[d_I][d_J] |= (0x01 << 4*d_p);
                              }
			    d_p++;
			    
			    if(d_p == 2){
			      d_p = 0;
			      ++d_J;
			      if(d_J == L/2){
				      d_J = 0;
				      ++d_I;
			      }
			      if(d_I < L && d_J < L/2)
				      Dir[d_I][d_J] = 0;
			   }
                           /*   
			   if(Score_Matrix[I+i][J+j] > sw_max){
				   A = I+i;
				   B = J+j;
				   sw_max = Score_Matrix[I+i][J+j];
			   }
			   */
                   } 
			
	       	}
            }
        } 
      }
              
        A = Score_Matrix[0][0];
	seq1_out[A] = 'Y';
	seq1_out[r] = 'C';
	seq1_out[size-r] = 'C';
	seq2_out[c] = 'D';
	seq2_out[size-c] = 'D';
	
	   
	 if(A >= B)
	      S_I = A;
         else
	      S_I = B;
      
	//seq1_out[A] = 'W';
        //seq2_out[B] = 'W';
	
   /*Traceback function*/
 /* 
     char SW_dir;
     uint8_t c1, c2;
     int count=0;
     uint8_t eo1, eo2;
     int prev_id = 0;
     int p_t;

     for(int n = L; n >=0; --n){
	if(prev_id != n/8){     
        s1_out[n/8] = 0;
        s2_out[n/8] = 0;
	}
	if(M[A][B]!=0 && n <= S_I){  
       		p_t = 1 - B%2;	
		if(A>=1 && B>=1){
		if(B%2 == 1)	
       		   SW_dir = (Dir[A-1][(B-1)/2] >> 4*p_t) & 0x0F;
		else
                   SW_dir = (Dir[A-1][B/2 - 1] >> 4*p_t) & 0x0F;
		
		eo1 = (seq1[A/8] >> 4*(A%8)) & 0x0F;
                eo2 = (seq2[B/8] >> 4*(B%8)) & 0x0F;

    		if(SW_dir == 'm'){
                        c1 = eo1;
    			c2 = eo2;
    			A = A-1;
    			B = B-1;
    		} else if(SW_dir == 'x'){
    		        c2 = 0x0D; 
    		   	c1 = eo1;
    		   	A = A-1;
    			}
    	       		else if(SW_dir == 'y'){
    	       	      		c1 = 0x0D;
    	       	      		c2 = eo2;
    	       	      		B = B-1;
    	            		}
		s1_out[n/8] |= (c1 << 4*(n%8));
	        s2_out[n/8] |= (c2 << 4*(n%8));
		}
		
       }
	else if(M[A][B] == 0  && n <=S_I){
		s1_out[n/8] |= (0x0E << ((n%8))*4);
	        s2_out[n/8] |= (0x0E << ((n%8))*4);
	     }else if(M[A][B] !=0 && n >S_I){
		s1_out[n/8] |= (0x0A << ((n%8))*4);
	        s2_out[n/8] |= (0x0A << ((n%8))*4);
	     }	
	     prev_id = n/8;
     
     }
     unpacking(s1_out, &seq1_out[sq_i]);
     unpacking(s2_out, &seq2_out[sq_i]);        
*/
    }
}


/*Main function*/
int main(int argc, char *argv[]){
    
    FILE *input1, *input2;
    FILE *output;
   /*Read in the two sequences to be aligned, one from refrence and another a query
    *short read, which are stored in a text file and store in arrays seq1[] and seq2[]
    */
    input1 = fopen("seq1_out.txt","rb");
	if (!input1) {
	  printf("Unable to open input file %s.\n", "seq1_out.txt");
	  fflush(stdout);
	  exit(-1);
	}	
	input2 = fopen("seq2_out.txt","rb");
	if (!input2) {
	  printf("Unable to open input file %s.\n", "seq2_out.txt");
	  fflush(stdout);
	  exit(-1);
	}

    output = fopen("align_out.txt","wb");
    
    char *seq1, *seq2;
    char *seq1_out, *seq2_out;
    char line[] = "Output seq 1:";
    char line1[] = "Output seq 2:";
    char head[] = "Sequence pair";
    int l_size = strlen(line);
    size_t  s_size = no_seq * (L+1) * sizeof(char) ;
    cudaEvent_t start, stop;
    float milliseconds;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
   

    /*Dynamic memory allocation at Host*/
    seq1 = (char*) malloc(s_size);
    if (seq1 == NULL) fprintf(stderr, "Bad malloc on seq1\n");
    seq2 = (char*) malloc(s_size);
    if (seq2 == NULL) fprintf(stderr, "Bad malloc on seq2\n");
    seq1_out = (char*) malloc(s_size);
    if (seq1_out == NULL) fprintf(stderr, "Bad malloc on seq1_out\n");
    seq2_out = (char*) malloc(s_size);
    if (seq2_out == NULL) fprintf(stderr, "Bad malloc on seq2_out\n");
 /*  
   int r_L;         //Reduced (L+1) is used as 8 sequence elements are stored in 1 entry, 4 bits each
   if((L+1)%8 != 0)
     r_L = ((L+1)/8 + 1);
   else
     r_L = (L+1)/8;

   size_t s_gpu = no_seq * r_L * sizeof(uint32_t);
  */ 
    /*Allocate memory in Device*/
    char *seq1_d;
    cudaMalloc(&seq1_d, s_size);
    char *seq2_d;
    cudaMalloc(&seq2_d, s_size);
    char *seq1_out_d;
    cudaMalloc(&seq1_out_d, s_size);
    char *seq2_out_d;
    cudaMalloc(&seq2_out_d, s_size);

    /* Load data from textfile */
    seq1[0] = '-';
    seq2[0] = '-';
    fread(&seq1[1], sizeof(char), ((L+1)*(no_seq-1)+ L), input1);
    fread(&seq2[1], sizeof(char), ((L+1)*(no_seq-1)+ L), input2);
     
   // printf("First char of seq1:%c, seq2:%c, last char of seq1:%c, seq2:%c\n", seq1[1], seq2[1], seq1[L], seq2[L]);

    fclose(input1);
    fclose(input2);
    fflush(stdout);
    
    cudaEventRecord(start);
    /*Copy data from Host to Device*/
    cudaMemcpy(seq1_d, seq1, s_size, cudaMemcpyHostToDevice);
    cudaMemcpy(seq2_d, seq2, s_size, cudaMemcpyHostToDevice);
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    printf("Time in H to D:%4.4f\n", milliseconds);

    
    /*Perform alignment at Device*/
    read_align<<<1,(no_seq)>>>(seq1_d, seq2_d, seq1_out_d, seq2_out_d);
  
    cudaDeviceSynchronize();
   
    cudaEventRecord(start);
    /*Copy output data from Device to Host*/
    cudaMemcpy(seq1_out, seq1_out_d, s_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(seq2_out, seq2_out_d, s_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(seq1, seq1_d, s_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(seq2, seq2_d, s_size, cudaMemcpyDeviceToHost);
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    printf("Time in D to H:%4.4f\n", milliseconds);

    /* Write result to file */
    for(int m=0; m < no_seq; m++){
	fwrite(head, sizeof(char), strlen(head), output);
        fprintf(output, "%d\n", m);	
        fwrite(line, sizeof(char), strlen(line), output);
        //fwrite(&seq1[m*(L+1)], sizeof(char), L+1, output);
        //fprintf(output,"\n");
        fwrite(&seq1_out[m*(L+1)], sizeof(char), L+1, output);
        fprintf(output,"\n");
        fwrite(line1, sizeof(char), strlen(line1), output);
        //fwrite(&seq2[m*(L+1)], sizeof(char), L+1, output);
        //fprintf(output, "\n");
        fwrite(&seq2_out[m*(L+1)], sizeof(char), L+1, output);
        if(m != no_seq-1)
          fprintf(output,"\n");
    }

	fclose(output);

	printf("Output complete.\n");
	fflush(stdout); 

	cudaError_t err= cudaGetLastError();
	if( err != cudaSuccess)
		printf("Error:%s\n",cudaGetErrorString(err));

    /*Free Device memory*/
    cudaFree(seq1_d);
    cudaFree(seq2_d);
    cudaFree(seq1_out_d);
    cudaFree(seq2_out_d);

    /*Free Host memory*/
    free(seq1);
    free(seq2);
    free(seq1_out);
    free(seq2_out);
}
