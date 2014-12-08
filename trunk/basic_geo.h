/** @file
 * Contains geometrical primitive operations.
 */
#ifndef BASIC_GEO_H
#define BASIC_GEO_H
#include <vector>

using std::vector;

namespace mesh_generation {
/**
 * Double precision error.
 */
extern const double kError;

/**
* Double comparison.
*/
int Cmp(double a, double b);

/**
* Returns the square of the value.
*/
long double Sqr(long double a);

struct Pixel;
/**
 * Euclidean 2D Point.
 */
struct Point {
	Point();
	Point(double x, double y);
	Point operator +(Point a) const;
	Point operator -(Point a) const;
	Point operator *(double a) const;
	Point operator /(double a) const;
	double operator *(Point a) const;

  /**
   * Returns 2D cross product.
   */
	double operator ^ (Point a) const;
	bool operator > (Point a) const;
	bool operator < (Point a) const;
	bool operator == (Point a) const;
  bool operator != (Point a) const;

  /**
   * Returns a ortogonal vector.
   */
	Point Ort() const;

  /**
   * @return the vector with the same direction, but modulus equal to one.
   */
  Point Unit() const;

  /**
   * Returns the distance to the Point a.
   */
	double Dist(Point a) const;

  /**
   * Transforms to pixel coordinates.
   */
  Pixel ToPixel(int image_height) const;

  /**
   * Coordinates.
   */
	double x;
  double y;
};

/**
 * Stores the coordinates of a pixel, the coordinates of pixels are bounded to
 * the image space [0, image_height-1] x [0, image_width-1].
 */
struct Pixel {
  Pixel();
  Pixel(int row, int col);

  /**
   * Transforms into cartesian coordinates.
   */
  Point ToPoint(int image_height) const;

  /**
   * Allows std::set insertion.
   */
  bool operator < (const Pixel& a) const;
  bool operator == (const Pixel& a) const;
  /**
   * Image coordinates.
   */
  int row;
  int col;
};

/**
 * Finds the circumcircle of thre tree Points a, b and c.
 *
 * @param center Returns the center of the circumcircle.
 * @param r Returns the radius of the circumcircle.
 */
void Circumcircle(Point a, Point b, Point c, Point* center , double* r);

/**
 * Checks if a point is strictly inside a circumcircle.
 *
 * @param d The point to be checked.
 * @param a,b,c The points of the triangle in counter-clockwise order.
 */
bool IsInsideCircumcircle(const Point& d, const Point& a, const Point& b,
                          const Point& c);

/**
 * Calculates the mean point of the polygon given in counter-clockwise order.
 */
Point MeanPoint(const vector<Point>& points);

/**
 * Calculates the centroid of the polygon given in counter-clockwise order.
 */
Point Centroid(const vector<Point>& points);

/**
 * Computes the barycentric coordinates (u, v, w) for point p with respect to
 * triangle (a, b, c).
 */
void Barycentric(const Point p, const Point a, const Point b, const Point c,
                 double* u, double* v, double* w);

}  // namespace mesh_generation
#endif
