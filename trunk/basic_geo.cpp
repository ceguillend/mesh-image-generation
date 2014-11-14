#include "basic_geo.h"
#include <math.h>
#include <assert.h>

namespace mesh_generation {

const double kError = 1e-7;

int Cmp(double a, double b) {
	if(fabs(a - b) < kError) return 0;
	if(a < b) return -1;
	return 1;
}

double Sqr(double a) {
	return a * a;
}

Point::Point(): x(0), y(0) {}

Point::Point(double x, double y): x(x), y(y) {}

double Point::Dist(Point a) const {
	return sqrt(Sqr(x - a.x) + Sqr(y - a.y));
}

double Point::operator^(Point a) const {
	return x * a.y - y * a.x;
}

Point Point::operator +(Point a) const {
	return Point(x+a.x, y+a.y);
}

Point Point::operator -(Point a) const {
	return Point(x-a.x, y-a.y);
}

double Point::operator* (Point a) const {
	return x * a.x + y * a.y;
}

Point Point::operator* (double a) const {
	return Point(x*a, y*a);
}

Point Point::operator/ (double a) const {
	return Point(x/a, y/a);
}

Point Point::Ort() const {
	return Point(-y, x);
}

bool Point::operator < (Point a) const {
  if (Cmp(x, a.x) != 0) return Cmp(x, a.x) < 0;
  return Cmp(y, a.y) < 0;
}

bool Point::operator > (Point a) const {
  if (Cmp(x, a.x) != 0) return Cmp(x, a.x) > 0;
  return Cmp(y, a.y) > 0;
}

bool Point::operator ==(Point a) const {
	return Cmp(x, a.x) == 0 && Cmp(y, a.y) == 0;
}

bool Point::operator != (Point a) const {
  return Cmp(x, a.x) != 0 || Cmp(y, a.y) != 0;
}


Pixel Point::ToPixel(int image_height) const{
	return Pixel(round(image_height - 1 - y), round(x));
}

Pixel::Pixel(): row(0), col(0) {}

Pixel::Pixel(int row, int col): row(row), col(col) {}

Point Pixel::ToPoint(int image_height) const {
  return Point(col, image_height - 1 - row);
}

bool Pixel::operator < (const Pixel& a) const {
  if (row != a.row) return row < a.row;
  return col < a.col;
}

bool Pixel::operator == (const Pixel& a) const {
  return row == a.row && col == a.col;
}

void Circumcircle(Point a, Point b, Point c, Point* center, double* r) {
	Point p = (a+b)/2, q = (b+c)/2;
	Point v = (b-a), u = (c-b);
  // Non co-linear points.
  assert(Cmp(v.Ort() * u, 0) != 0);
	double k = ((q-p)*u)/(v.Ort()*u);

	*center  = p + v.Ort()*k;
	*r = center->Dist(a);
}

bool IsInsideCircumcircle(const Point& d, const Point& a, const Point& b,
                          const Point& c)
{
	double A[3][3];
	A[0][0] = a.x-d.x; A[0][1] = a.y-d.y; A[0][2] =
      Sqr(a.x) - Sqr(d.x) + Sqr(a.y) - Sqr(d.y);

	A[1][0] = b.x-d.x; A[1][1] = b.y-d.y; A[1][2] =
      Sqr(b.x) - Sqr(d.x) + Sqr(b.y) - Sqr(d.y);

	A[2][0] = c.x-d.x; A[2][1] = c.y-d.y; A[2][2] =
      Sqr(c.x) - Sqr(d.x) + Sqr(c.y) - Sqr(d.y);

	double det = +A[0][0] * (A[1][1] * A[2][2] - A[2][1] * A[1][2])
	              -A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
		            +A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);

	return Cmp(det, 0) < 0;
}

Point Centroid(const vector<Point>& points) {
  Point ans;
  for (const Point& pt: points) {
    ans = ans + pt;
  }
  return ans/points.size();
}

Pixel Centroid(const vector<Pixel>& pixels) {
  Pixel ans;
  for (const Pixel pixel: pixels) {
    ans = Pixel(ans.row + pixel.row, ans.col + pixel.col);
  }
  ans.row /= pixels.size();
  ans.col /= pixels.size();
  return ans;
}
}  // namespace mesh_generation
