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

__device__ void init_DP(int M[][L+1], int X[][L+1], int Y[][L+1]){
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

__device__ void compute_DP(sw_entry SW_i_j, int seq1_i, int seq2_i, char *seq1, char *seq2, int M[][L+1], int X[][L+1], int Y[][L+1]){
    int M_max =0, X_max, Y_max;
    int match_score;
    //printf("BEFORE\n");
    //printf("Index i:%d, j:%d, seq1:%c, seq2:%c, score:%d, dir:%d\n",seq1_i, seq2_i, seq1[seq1_i], seq2[seq2_i], SW_i_j.value, SW_i_j.direction);
   if(seq1[seq1_i] == seq2[seq2_i])
   	match_score = match;
   else
   	match_score = mismatch;

	if(M[seq1_i-1][seq2_i-1] + match_score> M_max)
		M_max = M[seq1_i-1][seq2_i-1] + match_score;
	if(X[seq1_i-1][seq2_i-1] + match_score> M_max)
		M_max = X[seq1_i-1][seq2_i-1] + match_score;
	if(Y[seq1_i-1][seq2_i-1] + match_score> M_max)
		M_max = Y[seq1_i-1][seq2_i-1] + match_score;

	M[seq1_i][seq2_i] =  M_max;

    Y_max = gap_extn + Y[seq1_i][seq2_i-1];
    if(penalty + M[seq1_i][seq2_i-1] > Y_max)
    	Y_max = M[seq1_i][seq2_i-1] + penalty;

    Y[seq1_i][seq2_i] = Y_max;

    X_max = gap_extn + X[seq1_i-1][seq2_i];
    if(penalty + M[seq1_i-1][seq2_i] > X_max)
    	X_max = M[seq1_i-1][seq2_i] + penalty;

    X[seq1_i][seq2_i] = X_max;

    SW_i_j.value = M_max;
    SW_i_j.direction = m;
    if(SW_i_j.value < X_max){
    	SW_i_j.value = X_max;
    	SW_i_j.direction = x;
    }
    if(SW_i_j.value < Y_max){
    	SW_i_j.value = Y_max;
    	SW_i_j.direction = y;
    }
    //printf("AFTER\n");
    //printf("Index i:%d, j:%d, seq1:%c, seq2:%c, score:%d, dir:%d\n",seq1_i, seq2_i, seq1[seq1_i], seq2[seq2_i], SW_i_j.value, SW_i_j.direction);


}

//__device__ sw_entry sw_max;
//__device__ int idx_i, idx_j;
__device__ void traceback(sw_entry SW[][L+1], int M[][L+1], char *seq1, char *seq2, char *seq1_out, char *seq2_out){
	sw_entry sw_max;
	int idx_i, idx_j;

	sw_max = SW[0][0];
	for(int i=0; i < L+1; i++){
		for(int j=0; j < L+1; j++){
			if(SW[i][j].value > sw_max.value){
				sw_max = SW[i][j];
				idx_i = i;
				idx_j = j;
			}
		}
          }
    //printf("Highest score index i:%d, j:%d and score:%d\n",idx_i, idx_j, sw_max.value);
    int I = idx_i, J = idx_j;
    int s_idx;
    if(idx_i > idx_j)
    	s_idx = idx_i;
    else
    	s_idx = idx_j;
    seq1_out[s_idx+1] ='\0';
    seq2_out[s_idx+1] ='\0';

    while(M[I][J]){
        //printf("**Index I:%d, J:%d, s_idx:%d char in seq1:%c, seq2:%c\n", I, J, s_idx, seq1[I], seq2[J]);
    	if(SW[I][J].direction == m){
                seq1_out[s_idx] = seq1[I];
    		seq2_out[s_idx] = seq2[J];
    		I = I-1;
    		J = J-1;
    	} else if(SW[I][J].direction == x){
    		     seq2_out[s_idx] = '-';
    		     seq1_out[s_idx] = seq1[I];
    		     I = I-1;
    		   }
    	       else {
    	       	 seq1_out[s_idx] = '-';
    	       	 seq2_out[s_idx] = seq2[J];
    	       	 J = J-1;
    	       }
      //printf("Score of M: %d\n", M[I][J]);
      //printf("Index I:%d, J:%d, char in seq1_out:%c, seq2_out:%c\n", I, J, seq1_out[s_idx], seq2_out[s_idx]);
      --s_idx;
    }

    while(s_idx > 0){
    seq1_out[s_idx] = '*';
    seq2_out[s_idx] = '*';
    --s_idx;
    }

}

__global__ void read_align(char *seq1, char *seq2, char *seq1_out, char *seq2_out){
    
   int seq_i;
   sw_entry Score_Matrix[L+1][L+1];
   int M[L+1][L+1], X[L+1][L+1], Y[L+1][L+1];  //DP matrices

   int index = blockIdx.x * blockDim.x +threadIdx.x;

   if(index < no_seq)
   {   
        seq_i = index * (L+1);
        seq1[seq_i] = '-';
        seq2[seq_i] = '-';
        seq1_out[seq_i] = '$';
        seq2_out[seq_i] = '$';
        /*Start scoring*/
        
        init_DP(M, X, Y);
      
            Score_Matrix[0][0].value = 0;
            for(int j=1; j<L+1; j++){
              Score_Matrix[0][j].value = 0;
            }
            for(int i=1; i<L+1; i++){
              Score_Matrix[i][0].value = 0;
            }

        for(int i=1; i<L+1; i++){
            for(int j=1; j<L+1; j++){
                compute_DP(Score_Matrix[i][j], i,j, &seq1[seq_i], &seq2[seq_i], M, X, Y);
            }
        }

        traceback(Score_Matrix, M, &seq1[seq_i], &seq2[seq_i], &seq1_out[seq_i], &seq2_out[seq_i]);
              

    }
}


/*Main function*/
int main(int argc, char *argv[]){
    
    FILE *input1, *input2;
    FILE *output;
   /*Read in the two sequences to be aligned, one from refrence and another a query
    *short read, which are stored in a text file and store in arrays seq1[] and seq2[]
    */
    //sprintf(buff1,argv[1]);
    //sprintf(buff2,argv[2]);
    input1 = fopen("seq1_out.txt","rb");//input1 = fopen(argv[1],"rb");
	if (!input1) {
	  printf("Unable to open input file %s.\n", "seq1_out.txt");//argv[1]);
	  fflush(stdout);
	  exit(-1);
	}	
	input2 = fopen("seq2_out.txt","rb");//input2 = fopen(argv[2],"rb");
	if (!input2) {
	  printf("Unable to open input file %s.\n", "seq2_out.txt");//argv[2]);
	  fflush(stdout);
	  exit(-1);
	}

    output = fopen("align_out.txt","wb");
    
    char *seq1, *seq2;
    char *seq1_out, *seq2_out;
    char line[] = "Output seq 1:";
    char line1[] = "Output seq 2:";
    int l_size = strlen(line);
    size_t  s_size = no_seq * (L+1) * sizeof(char) ;
   

    /*Dynamic memory allocation at Host*/
    seq1 = (char*) malloc(s_size);
    if (seq1 == NULL) fprintf(stderr, "Bad malloc on seq1\n");
    seq2 = (char*) malloc(s_size);
    if (seq2 == NULL) fprintf(stderr, "Bad malloc on seq2\n");
    seq1_out = (char*) malloc(s_size);
    if (seq1_out == NULL) fprintf(stderr, "Bad malloc on seq1_out\n");
    seq2_out = (char*) malloc(s_size);
    if (seq2_out == NULL) fprintf(stderr, "Bad malloc on seq2_out\n");
   
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
     
    printf("First char of seq1:%c, seq2:%c, last char of seq1:%c, seq2:%c\n", seq1[1], seq2[1], seq1[L], seq2[L]);

    fclose(input1);
    fclose(input2);
    fflush(stdout);
    
    printf("Strlen of seq1:%d, seq2:%d\n", strlen(seq1), strlen(seq2));
    /*Copy data from Host to Device*/
    cudaMemcpy(seq1_d, seq1, s_size, cudaMemcpyHostToDevice);
    cudaMemcpy(seq2_d, seq2, s_size, cudaMemcpyHostToDevice);
   
    /*Perform alignment at Device*/
    read_align<<<1,1>>>(seq1_d, seq2_d, seq1_out_d, seq2_out_d);
  
    cudaDeviceSynchronize();
   
    /*Copy output data from Device to Host*/
    cudaMemcpy(seq1_out, seq1_out_d, s_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(seq2_out, seq2_out_d, s_size, cudaMemcpyDeviceToHost);
    printf("Strlen of seq1_out:%d, seq2_out:%d\n",strlen(seq1_out), strlen(seq2_out));
    /* Write result to file */
    for(int m=0; m < no_seq; m++){
        fwrite(line, sizeof(char), strlen(line), output);
        fwrite(&seq1[m*(L+1)], sizeof(char), L+1, output);
        fprintf(output,"\n");
        fwrite(&seq1_out[m*(L+1)], sizeof(char), L+1, output);
        fprintf(output,"\n");
        fwrite(line1, sizeof(char), strlen(line1), output);
        fwrite(&seq2[m*(L+1)], sizeof(char), L+1, output);
        fprintf(output, "\n");
        fwrite(&seq2_out[m*(L+1)], sizeof(char), L+1, output);
        if(m != no_seq-1)
          fprintf(output,"\n");
    }

	fclose(output);

	printf("Output complete.\n");
	fflush(stdout);

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

    return 0;
}
