#include"basic_geo.h"
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

/** Pixel coordinates to euclidean coordinates
*
*@Param px Pixel coordinates
$@Param h Height of the image
*@Return Euclidean coordinates
*/
point px_pt(pii px, int H) 
{
	return point(px.F, H-1-px.S);

}

/** Euclidean coordinates to pixel coordinates
*
*@Param pt Euclidean coordinates
$@Param h Height of the image
*@Return Pixel coordinates
*/
pii pt_px(point pt, int H)
{
	return pii(round(pt.x), round(H-1-pt.y));
}
