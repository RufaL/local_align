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

#define L 100
#define no_seq 10

/*Considering 300bp constitute a short read*/
int M[300][300], X[300][300], Y[300][300];  //DP matrices
/*const int for penalty*/
const int penalty = gap_open + gap_extn;    

__global__ int match_score(int i, int j, char *seq1, char *seq2){
	if(seq1[i] == seq2[j])
		return match;
	else
		return mismatch;
}

__global__ void init_DP(int seq1_len, int seq2_len){
	M[0][0] = 0;
	X[0][0] = -1000;
	Y[0][0] = -1000;
	for(int i=1; i <seq1_len; i++){
		M[i][0] = 0;
		X[i][0] = -1000;   //Just a large negative number
		Y[i][0] = -1000;
	}

	for(int j=1; j< seq2_len; j++){
		M[0][j] = 0;
		X[0][j] = -1000;   //Just a large negative number
		Y[0][j] = -1000;
	}
}

__global__ sw_entry compute_DP(int seq1_i, int seq2_i, char *seq1, char *seq2){
    int M_max =0, X_max, Y_max;
    sw_entry SW_i_j;
    //printf("BEFORE\n");
    //printf("Index i:%d, j:%d, seq1:%c, seq2:%c, score:%d, dir:%d\n",seq1_i, seq2_i, seq1[seq1_i], seq2[seq2_i], SW_i_j.value, SW_i_j.direction);

	if(M[seq1_i-1][seq2_i-1] + match_score(seq1_i, seq2_i, seq1, seq2)> M_max)
		M_max = M[seq1_i-1][seq2_i-1] + match_score(seq1_i, seq2_i, seq1, seq2);
	if(X[seq1_i-1][seq2_i-1] + match_score(seq1_i, seq2_i, seq1, seq2)> M_max)
		M_max = X[seq1_i-1][seq2_i-1] + match_score(seq1_i, seq2_i, seq1, seq2);
	if(Y[seq1_i-1][seq2_i-1] + match_score(seq1_i, seq2_i, seq1, seq2)> M_max)
		M_max = Y[seq1_i-1][seq2_i-1] + match_score(seq1_i, seq2_i, seq1, seq2);

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

    return SW_i_j;

}

__global__ void traceback(sw_entry *SW, int seq1_len, int seq2_len, char *seq1, char *seq2, char *seq1_out, char *seq2_out){
	sw_entry sw_max;
	int idx_i, idx_j;

	sw_max = SW[0][0];
	for(int i=0; i < seq1_len+1; i++){
		for(int j=0; j < seq2_len+1; j++){
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

    while(s_idx >= 0){
    seq1_out[s_idx] = '*';
    seq2_out[s_idx] = '*';
    --s_idx;
    }

}

__global__ void read_align(char *seq1, char *seq2){
    
   int seq_i;
   sw_entry Score_Matrix[L+1][L+1];

   if(threadIdx.x < no_seq)
   {
        seq_i = threadIdx.x * (L+1);
        seq1[seq_i] = '-';
        seq2[seq_i] = '-';
        seq1_out[seq_i] = '$';
        seq2_out[seq_i] = '$';
        /*Start scoring*/
        init_DP(L+1, L+1);
      
            Score_Matrix[0][0].value = 0;
            for(int j=1; j<L+1; j++){
              Score_Matrix[0][j].value = 0;
            }
            for(int i=1; i<L+1; i++){
              Score_Matrix[i][0].value = 0;
            }

        for(int i=1; i<L+1; i++){
            for(int j=1; j<L+1; j++){
                Score_Matrix[i][j] = compute_DP(i,j, &seq1[seq_i], &seq2[seq_i]);
            }
        }

        traceback(Score_Matrix, L, L, &seq1[seq_i], &seq2[seq_i], &seq1_out[seq_out_i], &seq2_out[seq_out_i]);
        

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

    /*Dynamic memory allocation*/
    seq1 = (char*) malloc(no_seq * (L+1) * sizeof(char));
    if (seq1 == NULL) fprintf(stderr, "Bad malloc on seq1\n");
    seq2 = (char*) malloc(no_seq * (L+1) * sizeof(char));
    if (seq2 == NULL) fprintf(stderr, "Bad malloc on seq2\n");
    seq1_out = (char*) malloc(no_seq * (L+1) * sizeof(char));
    if (seq1_out == NULL) fprintf(stderr, "Bad malloc on seq1_out\n");
    seq2_out = (char*) malloc(no_seq * (L+1) * sizeof(char));
    if (seq2_out == NULL) fprintf(stderr, "Bad malloc on seq2_out\n");

    /* Load data from textfile */
    seq1[0] = '-';
    seq2[0] = '-';
    fread(&seq1[1], sizeof(char), (L+1)*no_seq, input1);
    fread(&seq2[1], sizeof(char), (L+1)*no_seq, input2);

    fclose(input1);
    fclose(input2);
    fflush(stdout);

    read_align<<<1,no_seq>>>(seq1, seq2);

    /* Write result to file */
    for(int m=0; m < no_seq; m++){
        fwrite(line, sizeof(char), strlen(line), output);
        fwrite(&seq1_out[m*(L+1)], sizeof(char), L+1, output);
        fprintf(output,"\n");
        fwrite(line1, sizeof(char), strlen(line1), output);
        fwrite(&seq2_out[m*(L+1)], sizeof(char), L+1, output);
    }

	fclose(output);

	printf("Output complete.\n");
	fflush(stdout);

    /*Cleanup*/
    free(seq1);
    free(seq2);
    free(seq1_out);
    free(seq2_out);

	return 0;
}
