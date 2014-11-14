/** @file
 * Methods to produce the thinning of binary image.
 */
#ifndef THINNING_H
#define THINNING_H

#include "opencv2/imgproc/imgproc.hpp"

namespace mesh_generation {
/** Thinning sub-iteration 1
*
*@param source Binary image to be transformed (CV_8UC1).
*@param output Thinned binary image (CV_8UC1).
*/
void ThinningSubiteration1(const cv::Mat& source, cv::Mat* output); 

/** Thinning sub-iteration 2
*
*@param source Binary image to be transformed (CV_8UC1).
*@param output Thinned binary image (CV_8UC1).
*/
void ThinningSubiteration2(const cv::Mat& source, cv::Mat* output);

/** Thins a binary image (CV_8UC1).
*
*@param input_image Binary image to be trasnformed (CV_8UC1).
*@param output_image Thinned binary image (CV_8UC1).
*/
void ThinningProcedure(const cv::Mat& input_image, cv::Mat* output_image);
}  // namespace mesh_generation
#endif

