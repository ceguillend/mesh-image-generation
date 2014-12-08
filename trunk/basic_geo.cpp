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

long double Sqr(long double a) {
  return a * a;
}

Point::Point(): x(0), y(0) {}

Point::Point(double x, double y): x(x), y(y) {}

double Point::Dist(Point a) const {
  return sqrt(Sqr(x - a.x) + Sqr(y - a.y));
}

Point Point::Unit() const {
  double mod = sqrt(Sqr(x) + Sqr(y));
  return Point(x / mod, y / mod);
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
  if (Cmp(y, a.y) != 0) return Cmp(y, a.y) < 0;
  return Cmp(x, a.x) < 0;
}

bool Point::operator > (Point a) const {
  if (Cmp(y, a.y) != 0) return Cmp(y, a.y) > 0;
  return Cmp(x, a.x) > 0;
}

bool Point::operator ==(Point a) const {
  return Cmp(x, a.x) == 0 && Cmp(y, a.y) == 0;
}

bool Point::operator != (Point a) const {
  return Cmp(x, a.x) != 0 || Cmp(y, a.y) != 0;
}


Pixel Point::ToPixel(int image_height) const {
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
                          const Point& c) {
  double A[3][3];
  A[0][0] = a.x-d.x;
  A[0][1] = a.y-d.y;
  A[0][2] =
    Sqr(a.x) - Sqr(d.x) + Sqr(a.y) - Sqr(d.y);

  A[1][0] = b.x-d.x;
  A[1][1] = b.y-d.y;
  A[1][2] =
    Sqr(b.x) - Sqr(d.x) + Sqr(b.y) - Sqr(d.y);

  A[2][0] = c.x-d.x;
  A[2][1] = c.y-d.y;
  A[2][2] =
    Sqr(c.x) - Sqr(d.x) + Sqr(c.y) - Sqr(d.y);

  double det = +A[0][0] * (A[1][1] * A[2][2] - A[2][1] * A[1][2])
               -A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
               +A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);

  return Cmp(det, 0) < 0;
}

Point MeanPoint(const vector<Point>& points) {
  Point mean(0, 0);
  for (const Point& point: points) {
    mean = mean + point;
  }
  return mean/points.size();
}

Point Centroid(const vector<Point>& points) {
  assert (points.size() > 2);
  long double cx = 0;
  long double cy = 0;
  long double area = 0;
  for (int i = 0; i < points.size(); ++i) {
    int j = (i + 1) % points.size();
    long double x0 = points[i].x;
    long double y0 = points[i].y;
    long double x1 = points[j].x;
    long double y1 = points[j].y;
    long double qarea = x0 * y1 - x1 * y0;
    area += qarea;
    cx += (x0 + x1) * qarea;
    cy += (y0 + y1) * qarea;
  }
  area /= 2;
  cx /= 6 * area;
  cy /= 6 * area;
  return Point(cx, cy);
}

void Barycentric(const Point p, const Point a, const Point b, const Point c,
                 double* u, double* v, double* w) {
    const Point v0 = b - a;
    const Point v1 = c - a;
    const Point v2 = p - a;
    const double d00 = v0 * v0;
    const double d01 = v0 * v1;
    const double d11 = v1 * v1;
    const double d20 = v2 * v0;
    const double d21 = v2 * v1;
    const double denom = d00 * d11 - d01 * d01;
    *v = (d11 * d20 - d01 * d21) / denom;
    *w = (d00 * d21 - d01 * d20) / denom;
    *u = 1.0 - *v - *w;
}

}  // namespace mesh_generation
