//functions for transforming a binary edge image into a set of segment lines
#ifndef CURVES_H
#define CURVES_H

#include "header.h"
#include "basic_geo.h"
#include "color_utils.h"

// Defines the maximum value for the colors to display the curves.
extern const int MAX_LUM;

/**Given the canny image of edges finds the pixel components and simple curves
* by 8-neighborhood of pixels
*
*@param img_edge Image containing a Canny image (CV_8U1).
*@param min_curve_sz The minimum length size needed to be considered a curve(in
                      pixels).
*@param curves  Returns the list of curves as a set of consecutive pixel
                coordinates.
*@param img_group Returns the image of the pixel components (every component of
                  the same color).
*@param img_curve Return the image of the curves found (ever curve of the same
                  color).
*/
void find_pixel_curves( const cv::Mat& img_edge, int min_curve_len,
                        std::vector< std::vector<std::pii> >& px_curve, 
			                  cv::Mat& img_group, cv::Mat& img_curve);

/** Simplify the curve of pixels into a set of line segments (douglas peucker
*   algorithm)
*
* @param px_curve The set of curves to be simplified
* @param sg_curve Return the set of segment points (adjacent elements are
*                 neighboors)
* @param max_d Maximum distance allowed from a point to its segment
* @param max_len Maximum length allowed for a segment 
*/
void simplify_curve(std::vector< std::vector<std::pii> >& px_curve,
                    int H, int W, double max_d, double max_sz, 
		                std::vector< std::vector<point> >& sg_curve,
                    cv::Mat& img_curve_pt);

#endif

