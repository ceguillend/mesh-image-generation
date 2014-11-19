/** @file
 * Functions for transforming a binary edge image into a set of line segments.
 */
#ifndef CURVES_H
#define CURVES_H

#include "basic_geo.h"
#include "color_utils.h"
#include "opencv2/imgproc/imgproc.hpp"
#include <vector>
#include <list>


namespace mesh_generation {
/**
 * Given the thinned canny image of pixel-edges finds the pixel components and
 * simple curves by 8-neighborhood of pixels. It is a precondition that the
 * Canny edge has to have edged of thickness equal to 1.
 *
 * @param canny_edge Image containing a Thinned Canny image (CV_8U1).
 * @param min_curve_length The minimum length considered for a curve, measured
 *          in pixels.
 * @param point_curves Returns the list of curves as a set of consecutive
 *          cartesian point coordinates.
 * @param curves_plot Returns the image of the curves found, by drawing them
 *          with different colors.
 */
void FindPixelCurves(const cv::Mat& canny_edge, int min_component_size,
                      int curve_length_threshold,
                      std::list< std::vector<Point> >* point_curves,
                      cv::Mat* curves_plot);
/**
 * Simplifies the curve of pixels into a set of line segments (Douglas Peucker
 * algorithm).
 *
 * @param point_curves The set of curves to be simplified.
 * @param image_length
 * @param image_width
 * @param max_point_distance Maximum distance allowed from a Point to its
 *          segment in the curve, measured in pixels.
 * @param max_segment_length Maximum length allowed for a segment, measured in
 *          pixels.
 * @param simplified_curves Returns the simplified curves.
 * @param simplified_curves_plot Returns the curves plotted on the image domain.
 */
void SimplifyPointCurves(const std::list<std::vector<Point> >& point_curves,
                        int image_height, int image_width,
                        int max_point_distance, int max_segment_length,
                        std::list< std::vector<Point> >* simplified_curves,
                        cv::Mat* simplified_curves_plot);

void
}  // namespace mesh_generation
#endif
