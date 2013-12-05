#ifndef THINNING_H
#define THINNING_H

#include"header.h"

// Zhang-Suen algorithm
//functions to caltulate thinning of binary image
void thinSubiteration1(cv::Mat & pSrc, cv::Mat & pDst); 
void thinSubiteration2(cv::Mat & pSrc, cv::Mat & pDst);
void thinning_procedure(cv::Mat & inputarray, cv::Mat & outputarray);
#endif

