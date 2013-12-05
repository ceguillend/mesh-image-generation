#ifndef BASIC_GEO_H
#define BASIC_GEO_H
#include"header.h"
#define eps 1e-6

//double utilities
double cmp(double a, double b);
double sqr(double a);

// euclidean 2D point
class point
{
	public:
	double x,y;
	//contructors
	point();
	point(double x, double y);
	//operators
	double dist(point a) const;
	double operator^(point a) const;
	point operator -(point a) const;
	bool operator >(point a) const;
	bool operator <(point a) const;
};


// checks if a point is inside a circumcircle
bool inside_circumcircle(const point& d, const point& a, const point& b, const point& c);

//pixel coordinates transformation
point px_pt(std::pii px, int H);
std::pii pt_px(point pt, int H);

#endif
