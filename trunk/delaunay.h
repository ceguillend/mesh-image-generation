#ifndef DELAUNAY_H_
#define DELAUNAY_H_

#include"basic_geo.h"

class delaunay
{
	private:

	vector<point> pt; //poitns
	/** triangles (points ids)
	* Triples of points ids (clockwise order)
	*/
	vector< vector<int> > triangle; //triangles (points ids)
	/** neighbor triangle (triangle ids)
	* defines the neighbor considering the point position 
	* as the beginning of the edge in clockwise order 
	*/
	vector< vector<int> > neighbor; 
	/**
	* Defines if the id_tr is a leaf or not
	*/				
	vector< int > is_leaf;
	/** Tree of triangles
	* Adjacency list of the triangle tree structure (DAG)
	*/
	vector< vector<int> > tree; 

	//checking a point containment
	bool inside_triangle(int id, const point& pt);
	int border_triangle(int id, const& point pt);

	//tree functions
	int search_point(const point& pt);
	
	public :
	//constructos
	delaunay();
	delaunay(int W, int H);
	//updating functions
	void add_point(const point& pt);
	//ploting functions	
	void plot_points() const;

};

#endif
