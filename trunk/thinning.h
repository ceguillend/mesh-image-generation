#ifndef THINNING_H
#define THINNING_H

#include"header.h"

// Zhang-Suen algorithm
//functions to caltulate thinning of binary image
void thinSubiteration1(Mat & pSrc, Mat & pDst); 
void thinSubiteration2(Mat & pSrc, Mat & pDst);
void thinning_procedure(Mat & inputarray, Mat & outputarray);
#endif

