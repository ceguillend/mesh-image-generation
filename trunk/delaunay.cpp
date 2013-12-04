#include"delaunay.h"

delaunay::delaunay()
{
}

/**
* Points will lie inside [x:W,y:H] rectable
* 
*/
delaunay::delaunay(int W, int H)
{
	// big enough triangle
	pt.push_back( point(2*W,2*H) );
	pt.push_back( point(2*W,-4*H) );
	pt.push_back( point(-4*W,2*H) );

	// insert that triangle
	triangle.push_back( vector<int>(3,0) );

	// init outside triangle
	tr[0][0]=0; // clockwise
	tr[0][1]=1;
	tr[0][2]=2;
	// no neighbors -1 = exterior
	neighbor.push_back( vector<int>(3,-1) );

	
}


/**
* Check if a point p is inside-or in the border of the triangle in tr[id]
*/
bool delaunay::inside_triangle(int id, const point& p)
{
	for(int i=1;i<4;++i)
	{
		point& a = pt[triangle[id][i-1]], &b = pt[triangle[id][i%3]], 
		if( cmp( (b-a)^(p-a),0 ) < 0)
			return 0;
	}
	return 1;
}

/** Tests if the points belongs to a edge
*@returns the strat point of the edge in clockwise order,
*	  -1 if it is not at the edge
*/
int delaunay::border_triangle(int id, const point& p)
{
	for(int i=1;i<4;++i)
	{
		point& a = pt[triangle[id][i-1]], &b = pt[triangle[id][i%3]], 
		if( cmp( (b-a)^(p-a),0 ) == 0)
			return i-1;
	}
	return -1
}
/** Look for the containing triangle in the triangle structure
*
*@param point to the searched
*@return Triangle id (note that this a leaf in the triangle tree)
*	 If it belongs to an edge it returns any of the neighbouring triangles
*/
int delaunay::search_triangle(const point& p)
{
	int id = 0; // root

	while( !is_leaf[id] )
	{
		int nx = -1;
		for(int i=0;i<tree[id].size();++i) // for each triangle
			if( inside_triangle(tree[id][i], p) )
			{
				nx = tree[id][i];
				break;
			}
		assert(nx!=-1);
		id = nx;
		
	}
	return id;
}


/** Updates the neighbor a given triangle to a new one 
* This is called when changes happens on the other side and the triangle
* needs to know that information
*
*@param tr_id trinagle in which the information is going to change
*@param old_neigh triangle id of the old neighbor
*@param new_neigh triangle id of the new neighbor
*/
void delaunay::update_neighbor(int tr_id, int old_neigh, int new_neig)
{
	for(int i=0;i<3;++i)
		if( neighbor[tr_id][i] == old_neigh )
		{
			neighbor[tr_id][i] = new_neigh;
			return;
		}
	
	assert(0);
}
/** Finds the sharing edge
*@param tr_id triangle to search on
*@param neigh_id neighbor trinagle
*/
int neighbor_edge(int tr_id, int neigh_id)
{
	for(int i=0;i<3;++i)
		if(neighbor[tr_id][i] == neigh_id)
			return i;
	assert(0); // it always expects a answer
}
/** Checks if the triangle at the given edge satisfies the delaunay property
* if not it flip the edges and legalize again
*
*@param tr_id triangle to be checked
*#param ed_id edge to be checked
*/
void delaunay::legalize_edge(int tr_id, int ed_id)
{
	int ntr_id = neighbor[tr_id][ed_id];
	int ned_id = neighbor_edge(ntr_id, tr_id);

	vector<int>& tr = triangle[tr_id], & ntr = triangle[ntr_id];

	if( inside_circumcircle( pt[ntr[(ned_id+2)%3]], pt[tr[0]], pt[tr[1]], pt[tr[2]] ) )
	{ // ilegal

	}
	// do nothing
}
/* Adds a point to the border of a trinagle
*
*/
void delaunay::add_point_border(int tr_id, int pt_id, int ed_id)
{
	assert( is_leaf[tr_id] );
	// neighbor information 
	int ntr_id = neighbor[tr_id][ed_id];
	assert( is_leaf[ntr_id] );
	
	int ned_id = neighbor_edge(ntr_id, tr_id);
	//
	vector<int> tr(2),ed(2);
	tr[0] = tr_id; tr[1] = ntr_id;
	ed[0] = tr_id; ed[2] = ned_id;

	// create two sons per triangle
	for(int i=0;i<2;++i) // split each triangle
	{
		is_leaf[tr[i]] = 0;
		vector<int>& son = tree[tr[i]];
		for(int j=0;j<2;++j)
		{
			son.push_back(triangle.size());
			triangle.push_back(vector<int>(3,-1));
			neighbor.push_back(vector<int>(3,-1));
			tree.push_back(vector<int>());
			is_leaf.push_back(1);
		}
	}
	//update the neighborhood and triangle points
	for(int i=0;i<2;++i)
	{
		vector<int>& son = tree[tr[i]];
		vector<int>& n_son = tree[tr[i^1]];

		vector<int>& q_tr = triangle[tr[i]];
		vector<int>& q_ng = neighbor[tr[i]]; 
		
		// first son
		triangle[son[0]][0] = q_tr[ed[i]];
		triangle[son[0]][1] = pt_id;
		triangle[son[0]][2] = q_tr[(ed[i]+2)%3];

		neighbor[son[0]][0] = n_son[1];
		neighbor[son[0]][1] = son[1];
		neighbor[son[0]][2] = q_ng[(ed[i]+2)%3]
	
		update_neighbor(neighbor[son[0]][2], tr[i], son[0]);
		
		//second son
		triangle[son[1]][0] = pt_id;
		triangle[son[1]][1] = q_tr[(ed[i]+1)%3];
		triangle[son[1]][2] = q_tr[(ed[i]+2)%3];

		neighbor[son[1]][0] = n_son[0];
		neighbor[son[1]][1] = q_ng[(ed[i]+1)%3]; 
		neighbor[son[1]][2] = son[0];

		update_neighbor(neighbor[son[1]][1], tr[i], son[1]);	
	}
	//legalize edges
}

/** Insert a point inside a trinagle
* 
*/
void deaunay::add_point_inside(int tr_id, int pt_id)
{
	assert( is_leaf[tr_id] );
	vector<int>& son = tree[tr_id];
	assert(son.emtpy()); // should be empty
	
	// creates 3 new triangles
	is_leaf[tr_id] = 0;
	for(int i=0;i<n;++i)
	{
		son.push_back( triangle.size() );
		triangle.push_back( vector<int> (3,-1) );
		neighbor.push_back( vector<int> (3,-1) );
		tree.push_back(vector<int>());
		is_leaf.push_back( 1 );
	}
	// set the right values for the triangles
	for(int i=0;i<3;++i)
	{
		triangle[son[i]][0] = pt_id;
		neighbor[son[i]][0] = son[(i-1+3)%3];
		
		triangle[son[i]][1] = triangle[tr_id][i];
		neighbor[son[i]][1] = neighbor[tr_id][i];

		triangle[son[i]][2] = triangle[tr_id][(i+1)%3];
		neighbor[son[i]][2] = son[(i+1)%3];		
		
		// update the neighbor's information
		update_neighbor(neighbor[son[i]][1], tr_id, son[i]);
	}

	// legalize edges

}
/**
*@param pt Point to be inserted to the delaunay triangulatoin
*/
void delaunay::add_point(const point& pt)
{
	// id of the triangle containing the point
	int tr_id = search_triangle( pt );
	
	int pt_id = pt.size();
	pt.push_back(pt); // adds to the set of points,
			  // note that it assumes that points are unique

	int ed_id; // edge of this triangle containing the point 
	if( (ed_id = border_triangle(tr_id, pt)) != -1 )
	{ // border
		add_point_border(tr_id, pt_id, ed_id);	
	}
	else
	{ // inside,
		add_point_inside(tr_id, pt_id);
	}
			
}

void delaunay::plot_points()
{

	FILE* file = fopen("data_p","w");

	double minx,maxx,miny,maxy;

	minx = miny = 1./0.;
	maxx = maxy = -1./0.;
	//write points
	for(int i=0;i<pt.size();++i)
	{
		// point coordinates, type , size
		fprintf(file, "%lf\t%lf\n", pt[i].x, pt[i].y);
		minx = min( minx, pt[i].x);
		maxx = max( maxx, pt[i].x);
		miny = min( miny, pt[i].y);
		maxy = max( maxy, pt[i].y);
	}
	fclose(file);

	// start the ploting process
	FILE* pf = fopen("graph.conf", "w");
	
	fprintf(pf, "set xtic auto\nset ytic auto\n");
	fprintf(pf, "set title \"Delaunay\" \n set xlabel \"X\"\n");
	fprintf(pf, "set ylabel \"Y\" \nset pointsize 1\n");
	double xs = (maxx-minx)*0.2, ys = (maxy-miny)*0.2;
	fprintf(pf, "set xr [%lf : %lf]\n", minx-xs, maxx+xs);
	fprintf(pf, "set yr [%lf : %lf]\n", miny-ys, maxy+ys);
	fprintf(pf, "plot \"data_p\" using 1:2 with points pt 2 title \"Points\" ");

	fclose(pf);
	
	int result = system("gnuplot -persist graph.conf");
}
