/*
 * Header file for SW algorithm functions: scoring matrix computation, 
 * DP matrix recurrence function implementation, rules for scoring like
 * score for match/mismatch and affine-gap penalty score, data structure of 
 * entries in SW scoring matrix
 */
#ifndef __SWALGIN_H__
#define __SWALGIN_H__

#include "stdio.h"
#include "stdlib.h"
#include "stdint.h"

#define L 100
#define no_seq 1

/*Enumeration of DP matrices*/
enum DP_dir{
 	m,            //Value 0 for diagonal direction :match/mismatch
 	x,            //Value 1 for x-direction : gap in sequence x
 	y             //Value 2 for y-direction : gap in sequence y
};


/* Scoring Matrix entry data structure*/
typedef struct {
	int value;    //Value of SW matrix entry
    DP_dir direction;  //Direction: M, X or Y which giving maximum score
}sw_entry;

const int match = 10, mismatch = -2;      //Match and Mismatch scores
const int gap_open = -15, gap_extn = -7; //Gap opening and extension penalty, values example
 

//void init_DP(int M[][L+1], int X[][L+1], int Y[][L+1]);

//void compute_DP(sw_entry SW_i_j, int seq1_i, int seq2_i, char *seq1, char *seq2, int M[][L+1], int X[][L+1], int Y[][L+1]);

//void traceback(sw_entry SW[][L+1], int M[][L+1], char *seq1, char *seq2, char *seq1_out, char *seq2_out);

#endif
