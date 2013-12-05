#ifndef CURVES_H
#define CURVES_H

#include"header.h"
#include"basic_geo.h"
//functions for transforming a binary edge image into a set of segment lines

void find_pixel_curves( const cv::Mat& img_edge, int min_curve_len, std::vector< std::vector<std::pii> >& px_curve, 
			cv::Mat& img_group, cv::Mat& img_curve );
//
void simplify_curve( std::vector< std::vector<std::pii> >& px_curve, int H, int W, double max_d, double max_sz, 
		     std::vector< std::vector<point> >& sg_curve, cv::Mat& img_curve_pt);

#endif

