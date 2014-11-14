#include "color_utils.h"
#include <random>

using namespace cv;

Vec3b hsv_to_rgb(Vec3d hsv) {
  double r, g, b;
  double h = hsv[0], s = hsv[1], v = hsv[2];
  int i = int(h * 6);
  double f = h * 6 - i;
  double p = v * (1 - s);
  double q = v * (1 - f * s);
  double t = v * (1 - (1 - f) * s);

  switch(i % 6){
      case 0: r = v, g = t, b = p; break;
      case 1: r = q, g = v, b = p; break;
      case 2: r = p, g = v, b = t; break;
      case 3: r = p, g = q, b = v; break;
      case 4: r = t, g = p, b = v; break;
      case 5: r = v, g = p, b = q; break;
  }
  return Vec3b(r * 255,  g * 255, b * 255);
}

cv::Vec3b rand_color() {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0, 1);
  return hsv_to_rgb(Vec3d(dis(gen), 0.5, 0.95));
}
