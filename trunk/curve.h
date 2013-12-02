#ifndef CURVES_H
#define CURVES_H

#include"header.h"
#include"basic_geo.h"
//functions for transforming a binary edge image into a set of segment lines

void find_pixel_curves( const Mat& img_edge, int min_curve_len, vector< vector<pii> >& px_curve, Mat& img_group, Mat& img_curve );
void simplify_curve( vector< vector<pii> >& px_curve, int H, int W, double max_d, double max_sz, vector< vector<point> >& sg_curve, Mat& img_curve_pt);

#endif

