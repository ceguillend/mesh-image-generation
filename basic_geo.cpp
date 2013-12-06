#include"basic_geo.h"

using namespace std;
/**
* Double comparison
*/
double cmp(double a, double b)
{
	if(fabs(a-b)<eps) return 0;
	if(a < b) return -1;
	return 1;
}
/**
* Square
*/
double sqr(double a)
{
	return a*a;
}

/**
* Point class
*/

point::point(){}
//
point::point(double x, double y):x(x),y(y){}
//
double point::dist(point a) const
{
	return sqrt(sqr(x-a.x)+sqr(y-a.y));
}
//
double point::operator^(point a) const
{
	return x*a.y-y*a.x;
}
point point::operator +(point a) const
{
	return point(x+a.x,y+a.y);
}
point point::operator -(point a) const
{
	return point(x-a.x,y-a.y);
}
double point::operator* (point a) const
{
	return x*a.x+y*a.y;
}
point point::operator* (double a) const
{
	return point(x*a, y*a);
}
point point::operator/ (double a) const
{
	return point(x/a, y/a);
}
point point::ort() const
{
	return point(-y, x);
}
bool point::operator <(point a) const
{
	return cmp(y,a.y)<0 || (cmp(y,a.y)==0 && cmp(x,a.x)<0);
}
bool point::operator >(point a) const
{
	return cmp(y,a.y)>0 || (cmp(y,a.y)==0 && cmp(x,a.x)>0);
}
bool point::operator ==(point a) const
{
	return cmp(x,a.x)==0 && cmp(y,a.y)==0;
}
/////


void circumcircle(point a, point b, point c, point& ce , double& r)
{
	point p = (a+b)/2, q = (b+c)/2;
	point v = (b-a), u = (c-b);
	double k = ((q-p)*u)/(v.ort()*u);
	
	ce  = p + v.ort()*k;
	r = ce.dist(a);
}

/** Check if the point is inside the circumcircle
*@param d is the point to be checked
*@param a,b,c are the points of the triangle in clockwise order
*/
bool inside_circumcircle(const point& d, const point& a, const point& b, const point& c)
{
	double A[3][3];

	A[0][0] = a.x-d.x; A[0][1] = a.y-d.y; A[0][2] = sqr(a.x)-sqr(d.x)+sqr(a.y)-sqr(d.y);
	A[1][0] = b.x-d.x; A[1][1] = b.y-d.y; A[1][2] = sqr(b.x)-sqr(d.x)+sqr(b.y)-sqr(d.y);
	A[2][0] = c.x-d.x; A[2][1] = c.y-d.y; A[2][2] = sqr(c.x)-sqr(d.x)+sqr(c.y)-sqr(d.y);

	double det = +A[0][0]*(A[1][1]*A[2][2]-A[2][1]*A[1][2])
	        -A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])
		+A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);

	return cmp(det,0) < 0;
}

/** Pixel coordinates to euclidean coordinates
*
*@Param px Pixel coordinates
$@Param h Height of the image
*@Return Euclidean coordinates
*/
point px_pt(pii px, int H) 
{
	return point(px.S, H-1-px.F);
		   // col,   row

}

/** Euclidean coordinates to pixel coordinates
*
*@Param pt Euclidean coordinates
$@Param h Height of the image
*@Return Pixel coordinates
*/
pii pt_px(point pt, int H)
{
	return pii( round(H-1-pt.y), round(pt.x));
		// row 			col
}
