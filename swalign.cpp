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

/*Considering 300bp constitute a short read*/
int M[300][300], X[300][300], Y[300][300];  //DP matrices

int match_score(int i, int j, char *seq1, char *seq2){
	if(seq1[i] == seq2[j])
		return match;
	else
		return mismatch;
}

void init_DP(int seq1_len, int seq2_len){
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

sw_entry compute_DP(int seq1_i, int seq2_i, char *seq1, char *seq2){
    int M_max =0, X_max, Y_max;
    int penalty = gap_open + gap_extn;
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

void traceback(sw_entry SW[301][301], int seq1_len, int seq2_len, char *seq1, char *seq2, char *seq1_out, char *seq2_out){
	sw_entry sw_max;
	int idx_i, idx_j;

	sw_max = SW[0][0];
	for(int i=0; i < seq1_len; i++){
		for(int j=0; j < seq2_len; j++){
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

#define L 100
#define no_seq 10

/*Main function*/
int main(int argc, char *argv[]){
    char seq1[301], seq2[301];
    char seq1_out[301], seq2_out[301];
    sw_entry Score_Matrix[301][301];
    int l1, l2;
    char buff1[128], buff2[128];
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

    char line[] = "Output seq 1:";
    char line1[] = "Output seq 2:";
    output = fopen("align_out.txt","wb");

    for (int k=0; k < no_seq; k++)
    {
	/* Load data from textfile */
        seq1[0] = '-';
        seq2[0] = '-';
	fread(&seq1[1], sizeof(char), L+1, input1);
	fread(&seq2[1], sizeof(char), L+1, input2); 
    
 
 	l1 = strlen(seq1)-1;
 	l2 = strlen(seq2)-1;
        printf("Size of seq1:%d\n", l1-1);
        printf("first char in seq1:%c\n", seq1[1]);
        printf("Size of seq2:%d\n", l2-1);
        printf("first char in seq2:%c\n", seq2[1]);
 	/*Start scoring*/
 	init_DP(l1, l2);
  
        Score_Matrix[0][0].value = 0;
        for(int j=1; j<L; j++){
          Score_Matrix[0][j].value = 0;
        }
        for(int i=1; i<L; i++){
          Score_Matrix[i][0].value = 0;
        }

 	for(int i=1; i<L; i++){
 		for(int j=1; j<L; j++){
 			Score_Matrix[i][j] = compute_DP(i,j, seq1, seq2);
 		}
 	}

    traceback(Score_Matrix, l1, l2, seq1, seq2, seq1_out, seq2_out);
    //printf("output seq1[1]:%c, seq2[1]:%c\n", seq1_out[1], seq2_out[1]);
    //printf("output length 1: %d, 2: %d\n", strlen(seq1_out), strlen(seq2_out));
    /* Write result to file */
      
        fwrite(line, sizeof(char), strlen(line), output);
	fwrite(seq1_out, sizeof(char), strlen(seq1_out), output);
        fprintf(output,"\n");
        fwrite(line1, sizeof(char), strlen(line1), output);
	fwrite(seq2_out, sizeof(char), strlen(seq2_out), output);
        if(k != L-1)
          fprintf(output,"\n");

    }
        fclose(input1);
        fclose(input2);
	fflush(stdout);

	fclose(output);

	printf("Output complete.\n");
	fflush(stdout);

	return 0;
}
