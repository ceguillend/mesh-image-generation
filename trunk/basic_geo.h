#ifndef BASIC_GEO_H
#define BASIC_GEO_H
#include"header.h"
#define eps 1e-6

//double utilities
double cmp(double a, double b);
double sqr(double a);

// euclidean 2D point
struct point
{
	double x,y;

	point(){}
	point(double x, double y):x(x),y(y){}
	double dist(point a)
	{
		return sqrt(sqr(x-a.x)+sqr(y-a.y));
	}
};

//pixel coordinates transformation
point px_pt(pii px, int H);
pii pt_px(point pt, int H);

#endif
