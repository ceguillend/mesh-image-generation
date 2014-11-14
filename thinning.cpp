#include"thinning.h"
#include <algorithm>
using namespace cv;

using cv::Mat;
using std::min;

namespace mesh_generation {

void ThinningSubiteration1(const cv::Mat & source, cv::Mat * output) {
  assert(source.type() == CV_8UC1);
  assert(output->type() == CV_8UC1);

  int rows = source.rows;
  int cols = source.cols;
  source.copyTo(*output);
  for(int i = 1; i < rows-1; i++) {
    for(int j = 1; j < cols-1; j++) {
      if(source.at<uchar>(i, j) == 255) {
        // Gets 8 neighbors.
        // Calculates C(p).
        int neighbor0 = source.at<uchar>(i-1, j-1) > 0;
        int neighbor1 = source.at<uchar>(i-1, j  ) > 0;
        int neighbor2 = source.at<uchar>(i-1, j+1) > 0;
        int neighbor3 = source.at<uchar>(i  , j+1) > 0;
        int neighbor4 = source.at<uchar>(i+1, j+1) > 0;
        int neighbor5 = source.at<uchar>(i+1, j  ) > 0;
        int neighbor6 = source.at<uchar>(i+1, j-1) > 0;
        int neighbor7 = source.at<uchar>(i  , j-1) > 0;

        int C = int(~neighbor1 & ( neighbor2 | neighbor3)) +
                int(~neighbor3 & ( neighbor4 | neighbor5)) +
                int(~neighbor5 & ( neighbor6 | neighbor7)) +
                int(~neighbor7 & ( neighbor0 | neighbor1));
        if(C == 1) {
          // Calculates N.
          int N1 = int(neighbor0 | neighbor1) +
                   int(neighbor2 | neighbor3) +
                   int(neighbor4 | neighbor5) +
                   int(neighbor6 | neighbor7);
          int N2 = int(neighbor1 | neighbor2) +
                   int(neighbor3 | neighbor4) +
                   int(neighbor5 | neighbor6) +
                   int(neighbor7 | neighbor0);
          int N = min(N1, N2);
          if ((N == 2) || (N == 3)) {
            /// Calculates criteria 3.
            int c3 = ( neighbor1 | neighbor2 | ~neighbor4) & neighbor3;
            if(c3 == 0) {
              output->at<uchar>(i, j) = 0;
            }
          }
        }
      }
    }
  }
}


void ThinningSubiteration2(const cv::Mat & source, cv::Mat * output) {
  assert(source.type() == CV_8UC1);
  assert(output->type() == CV_8UC1);

  int rows = source.rows;
  int cols = source.cols;
  source.copyTo(*output);

  for(int i = 1; i < rows-1; i++)  {
    for(int j = 1; j < cols-1; j++) {
      if (source.at<uchar>( i, j) == 255) {
        /// Gets 8 neighbors.
        /// Calculates C(p).
        int neighbor0 = source.at<uchar>(i-1, j-1) > 0;
        int neighbor1 = source.at<uchar>(i-1, j  ) > 0;
        int neighbor2 = source.at<uchar>(i-1, j+1) > 0;
        int neighbor3 = source.at<uchar>(i  , j+1) > 0;
        int neighbor4 = source.at<uchar>(i+1, j+1) > 0;
        int neighbor5 = source.at<uchar>(i+1, j  ) > 0;
        int neighbor6 = source.at<uchar>(i+1, j-1) > 0;
        int neighbor7 = source.at<uchar>(i  , j-1) > 0;
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
          int N = min(N1, N2);
          if((N == 2) || (N == 3)) {
            int E = (neighbor5 | neighbor6 | ~neighbor0) & neighbor7;
            if(E == 0) {
              output->at<uchar>(i, j) = 0;
            }
          }
        }
      }
    }
  }
}


namespace {
bool Equal(const cv::Mat& a, const cv::Mat& b) {
  assert(a.rows == b.rows && a.cols == b.cols);
  for (int i = 0; i < a.rows; ++i) {
    for (int j = 0; j < a.cols; ++j) {
      if (a.at<uchar>(i, j) != b.at<uchar>(i, j)) {
        return false;
      }
    }
  }
  return true;
}
}  // namespace

void ThinningProcedure(const cv::Mat& input_image, cv::Mat* output_image) {
  assert(input_image.type() == CV_8UC1);

  int rows = input_image.rows;
  int cols = input_image.cols;

  *output_image = Mat::zeros(rows, cols, CV_8UC1);

  // Enlarges input_image.
  Mat enlarged_source = Mat::zeros(rows + 2, cols + 2, CV_8UC1);

  for(int i = 0; i < rows; i++)
    for(int j = 0; j < cols; j++)
      enlarged_source.at<uchar>(i+1, j+1) = input_image.at<uchar>(i, j);

  Mat thinned_image_1 = Mat::zeros(rows + 2, cols + 2, CV_8UC1);
  Mat thinned_image_2 = Mat::zeros(rows + 2, cols + 2, CV_8UC1);

  int iterations = 0;
  while (true) {
    ThinningSubiteration1(enlarged_source, &thinned_image_1);
    ThinningSubiteration2(thinned_image_1, &thinned_image_2);
    ++iterations;
    if (Equal(enlarged_source, thinned_image_2)) {
      break;
    }
    thinned_image_2.copyTo(enlarged_source);
  }
  printf("Thinning iterations: %d\n", iterations); 

  // Writes the output.
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      output_image->at<uchar>(i, j) = enlarged_source.at<uchar>(i + 1, j + 1);
    }
  }
}
}  // namespace mesh_generation
