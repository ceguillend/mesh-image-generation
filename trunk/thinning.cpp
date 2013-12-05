#include"thinning.h"
using namespace cv;

/** Thining sub-iteration 1
*
*@param pSrc Binary image to be trasnformed (CV_8UC1)
*@param pDst Thinned binary image (CV_8UC1)
*/
void thinSubiteration1(Mat & pSrc, Mat & pDst) 
{
	assert(pSrc.type() == CV_8UC1);
	assert(pDst.type() == CV_8UC1);

        int rows = pSrc.rows;
        int cols = pSrc.cols;
        pSrc.copyTo(pDst);
        for(int i = 0; i < rows; i++) {
                for(int j = 0; j < cols; j++) {
                        if(pSrc.at<uchar>(i, j) == 255) { // marked 1
                                /// get 8 neighbors
                                /// calculate C(p)
                                int neighbor0 = pSrc.at<uchar>( i-1, j-1) > 0;
                                int neighbor1 = pSrc.at<uchar>( i-1, j  ) > 0;
                                int neighbor2 = pSrc.at<uchar>( i-1, j+1) > 0;
                                int neighbor3 = pSrc.at<uchar>( i  , j+1) > 0;
                                int neighbor4 = pSrc.at<uchar>( i+1, j+1) > 0;
                                int neighbor5 = pSrc.at<uchar>( i+1, j  ) > 0;
                                int neighbor6 = pSrc.at<uchar>( i+1, j-1) > 0;
                                int neighbor7 = pSrc.at<uchar>( i  , j-1) > 0;

                                int C = int(~neighbor1 & ( neighbor2 | neighbor3)) +
                                                 int(~neighbor3 & ( neighbor4 | neighbor5)) +
                                                 int(~neighbor5 & ( neighbor6 | neighbor7)) +
                                                 int(~neighbor7 & ( neighbor0 | neighbor1));
                                if(C == 1) {
                                        /// calculate N
                                        int N1 = int(neighbor0 | neighbor1) +
                                                         int(neighbor2 | neighbor3) +
                                                         int(neighbor4 | neighbor5) +
                                                         int(neighbor6 | neighbor7);
                                        int N2 = int(neighbor1 | neighbor2) +
                                                         int(neighbor3 | neighbor4) +
                                                         int(neighbor5 | neighbor6) +
                                                         int(neighbor7 | neighbor0);
                                        int N = min(N1,N2);
                                        if ((N == 2) || (N == 3)) {
                                                /// calculate criteria 3
                                                int c3 = ( neighbor1 | neighbor2 | ~neighbor4) & neighbor3;
                                                if(c3 == 0) {
                                                        pDst.at<uchar>( i, j) = 0;
                                                }
                                        }
                                }
                        }
                }
        }
}

/** Thining sub-iteration 2
*
*@param pSrc Binary image to be trasnformed (CV_8UC1)
*@param pDst Thinned binary image (CV_8UC1)
*/
void thinSubiteration2(Mat & pSrc, Mat & pDst) 
{
	assert(pSrc.type() == CV_8UC1);
	assert(pDst.type() == CV_8UC1);

        int rows = pSrc.rows;
        int cols = pSrc.cols;
        pSrc.copyTo( pDst);

        for(int i = 0; i < rows; i++) {
                for(int j = 0; j < cols; j++) {
                        if (pSrc.at<uchar>( i, j) == 255) { // marked 1
                                /// get 8 neighbors
                                /// calculate C(p)
                            int neighbor0 = pSrc.at<uchar>( i-1, j-1) > 0;
                            int neighbor1 = pSrc.at<uchar>( i-1, j  ) > 0;
                            int neighbor2 = pSrc.at<uchar>( i-1, j+1) > 0;
                            int neighbor3 = pSrc.at<uchar>( i  , j+1) > 0;
                            int neighbor4 = pSrc.at<uchar>( i+1, j+1) > 0;
                            int neighbor5 = pSrc.at<uchar>( i+1, j  ) > 0;
                            int neighbor6 = pSrc.at<uchar>( i+1, j-1) > 0;
                            int neighbor7 = pSrc.at<uchar>( i  , j-1) > 0;
                                int C = int(~neighbor1 & ( neighbor2 | neighbor3)) +
                                        int(~neighbor3 & ( neighbor4 | neighbor5)) +
                                        int(~neighbor5 & ( neighbor6 | neighbor7)) +
                                        int(~neighbor7 & ( neighbor0 | neighbor1));
                                if(C == 1) {
                                        /// calculate N
                                        int N1 = int(neighbor0 | neighbor1) +
                                                int(neighbor2 | neighbor3) +
                                                int(neighbor4 | neighbor5) +
                                                int(neighbor6 | neighbor7);
                                        int N2 = int(neighbor1 | neighbor2) +
                                                int(neighbor3 | neighbor4) +
                                                int(neighbor5 | neighbor6) +
                                                int(neighbor7 | neighbor0);
                                        int N = min(N1,N2);
                                        if((N == 2) || (N == 3)) {
                                                int E = (neighbor5 | neighbor6 | ~neighbor0) & neighbor7;
                                                if(E == 0) {
                                                        pDst.at<uchar>(i, j) = 0;
                                                }
                                        }
                                }
                        }
                }
        }
}

/** Thins a binary image (CV_8UC1)
*
*@param inputarray Binary image to be trasnformed (CV_8UC1)
*@param outputarray  Thinned binary image (CV_8UC1)
*/

void thinning_procedure(Mat & inputarray, Mat & outputarray)
{

	assert(inputarray.type() == CV_8UC1);

        int rows = inputarray.rows;
        int cols = inputarray.cols;

	outputarray = Mat::zeros(rows , cols , CV_8UC1);

        /// enlarge the inputarray 
        Mat p_enlarged_src = Mat(rows + 2, cols + 2, CV_8UC1);

        for(int i = 0; i < rows; i++) 
                for(int j = 0; j < cols; j++) 
                                p_enlarged_src.at<uchar>( i+1, j+1) = inputarray.at<uchar>(i, j);

        /// start to thin
        Mat p_thinMat1 = Mat::zeros(rows + 2, cols + 2, CV_8UC1);
        Mat p_thinMat2 = Mat::zeros(rows + 2, cols + 2, CV_8UC1);
        Mat p_cmp = Mat::zeros(rows + 2, cols + 2, CV_8UC1);

        while (1) {
                /// sub-iteration 1
                thinSubiteration1(p_enlarged_src, p_thinMat1);
                /// sub-iteration 2
                thinSubiteration2(p_thinMat1, p_thinMat2);
                /// compare
                compare(p_enlarged_src, p_thinMat2, p_cmp, CV_CMP_EQ);
		// check
                if(countNonZero(p_cmp) == (rows + 2) * (cols + 2))  //all the same
		{
                       break;
                }
                /// copy
                p_thinMat2.copyTo(p_enlarged_src);
        }

        // copy result
        for(int i = 0; i < rows; i++) 
                for(int j = 0; j < cols; j++) 
                        outputarray.at<uchar>(i, j) = p_enlarged_src.at<uchar>(i+1, j+1);
                
}
