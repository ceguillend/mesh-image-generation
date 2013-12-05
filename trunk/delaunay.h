#ifndef DELAUNAY_H_
#define DELAUNAY_H_

#include"basic_geo.h"

class delaunay
{
	private:

	std::vector<point> pt; //poitns
	/** triangles (points ids)
	* Triples of points ids (clockwise order)
	*/
	std::vector< std::vector<int> > triangle; //triangles (points ids)
	/** neighbor triangle (triangle ids)
	* defines the neighbor considering the point position 
	* as the beginning of the edge in clockwise order 
	*/
	std::vector< std::vector<int> > neighbor; 
	/** Tree of triangles
	* Adjacency list of the triangle tree structure (DAG)
	*/
	std::vector< std::vector<int> > tree; 
	/**
	* Number of triangles visited in all the queries
	*/
	int n_location_operations;

	//checking a point containment
	bool inside_triangle(int id, const point& pt) const;
	int border_triangle(int id, const point& pt) const;

	//tree functions
	int search_triangle(const point& pt);// const;

	// update functions
	int new_triangle();	
	void update_neighbor(int tr_id, int old_nei, int new_nei);	
	void legalize_edge(int tr_id, int ed_id);
	void add_point_border(int tr_id, int pt_id, int ed_id);
	void add_point_inside(int tr_id, int pt_id);
	// query functions
	int neighbor_edge(int tr_id, int neigh_id) const;

	public :
	//constructos
	delaunay();
	delaunay(int W, int H);
	//updating functions
	void add_point(point pt);
	// check correctness
	bool check() const;
	//ploting functions	
	void plot_points() const;
	void plot_triangulation() const;
	// statistics
	int size() const;
	double average_location_operations() const;
};

#endif
