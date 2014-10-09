// Defines functions for color generation
#ifndef COLOR_UTILS_H
#define COLOR_UTILS_H

#include "opencv2/imgproc/imgproc.hpp"

// transforms a hsv value into rgb format.
cv::Vec3b hsv_to_rgb(cv::Vec3d hsv);

// Generates a bright random color.
cv::Vec3b rand_color();

#endif // COLOR_UTILS_H
