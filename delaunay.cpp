#include"delaunay.h"

using namespace std;

delaunay::delaunay()
{
}

/** Allocates space for a new triangle
*@return The id of the allocated triangle
*/
int delaunay::new_triangle()
{
	int t = triangle.size();
	triangle.push_back(vector<int>(3,-1));
	neighbor.push_back(vector<int>(3,-1)); //exterior negihbors
	tree.push_back(vector<int>());
	return t;
}

/**
* Points will lie inside [x:0-W,y:0-H] rectable
* 
*/
delaunay::delaunay(int W, int H)
{
	// it assumed that (0,0) will be always part of the triangulation

	// big enough triangle
	pt.push_back( point(2*W,2*H) );
	pt.push_back( point(2*W,-4*H) );
	pt.push_back( point(-4*W,2*H) ); // part of the set points

	// insert that triangle
	int t = new_triangle();
	// init outside triangle
	for(int i=0;i<3;++i)
		triangle[t][i]=i; // clockwise

}


/**
* Check if a point p is inside-or in the border of the triangle in tr[id]
*/
bool delaunay::inside_triangle(int id, const point& p) const
{
	for(int i=1;i<4;++i)
	{
		const point& a = pt[triangle[id][i-1]],&b = pt[triangle[id][i%3]]; 
		if( cmp( (b-a)^(p-a),0 ) > 0)
			return 0;
	}
	return 1;
}

/** Tests if the points belongs to a edge
*@returns the strat point of the edge in clockwise order,
*	  -1 if it is not at the edge
*/
int delaunay::border_triangle(int id, const point& p) const
{
	for(int i=1;i<4;++i)
	{
		const point& a = pt[triangle[id][i-1]], &b = pt[triangle[id][i%3]]; 
		if( cmp( (b-a)^(p-a),0 ) == 0)
			return i-1;
	}
	return -1;
}
/** Look for the containing triangle in the triangle structure
*
*@param point to the searched
*@return Triangle id (note that this a leaf in the triangle tree)
*	 If it belongs to an edge it returns any of the neighbouring triangles
*/
int delaunay::search_triangle(const point& p) const
{
	int id = 0; // root

	while( !tree[id].empty() )
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
void delaunay::update_neighbor(int tr_id, int old_neigh, int new_neigh)
{
	if( tr_id == -1 )
		return;
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
int delaunay::neighbor_edge(int tr_id, int neigh_id) const
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
*#param ed_id edge to be checked( check the oposite triangle)
*/
void delaunay::legalize_edge(int tr_id, int ed_id)
{
	int ntr_id = neighbor[tr_id][ed_id];
	if( ntr_id == -1 ) // outside
		return ;
	int ned_id = neighbor_edge(ntr_id, tr_id);

	vector<int>& tri = triangle[tr_id], & ntri = triangle[ntr_id];

	if( inside_circumcircle( pt[ntri[(ned_id+2)%3]], pt[tri[0]], pt[tri[1]], pt[tri[2]] ) )
	{ // ilegal
		vector<int> tr(2,0), ed(2,0);
		tr[0] = tr_id;   ed[0] = ed_id;
		tr[1] = ntr_id; ed[1] = ned_id;
		//FLip edges
		vector<int> son(2);
		for(int i=0;i<2;++i)
			son[i] = new_triangle();
		// add the sons
		for(int i=0;i<2;++i)
			tree[tr[i]] = son;

		// udate adjecy and point values
		for(int i=0;i<2;++i)
		{
			triangle[son[i]][0] = triangle[tr[i]][(ed[i]+2)%3];
			triangle[son[i]][1] = triangle[tr[i]][ed[i]];
			triangle[son[i]][2] = triangle[tr[i^1]][(ed[i^1]+2)%3];

			neighbor[son[i]][0] = neighbor[tr[i]][(ed[i]+2)%3];
			neighbor[son[i]][1] = neighbor[tr[i^1]][(ed[i^1]+1)%3];
			neighbor[son[i]][2] = son[i^1];

			update_neighbor(neighbor[son[i]][0], tr[i], son[i]);
			update_neighbor(neighbor[son[i]][1], tr[i^1], son[i]);
		}
		//legalize
		//for(int i=0;i<2;++i)
		legalize_edge(son[0], 1);
		legalize_edge(son[1], 0);
	}
	// do nothing
}
/* Adds a point to the border of a trinagle
*
*/
void delaunay::add_point_border(int tr_id, int pt_id, int ed_id)
{
	assert( tree[tr_id].empty() ); // leaf
	// neighbor information 
	int ntr_id = neighbor[tr_id][ed_id];
	assert( tree[ntr_id].empty() ); // leaf
	
	int ned_id = neighbor_edge(ntr_id, tr_id);
	//
	vector<int> tr(2),ed(2);
	tr[0] = tr_id; tr[1] = ntr_id;
	ed[0] = ed_id; ed[1] = ned_id;

	// create two sons per triangle
	for(int i=0;i<2;++i) // split each triangle
	{
	//	vector<int>& son = tree[tr[i]];
		for(int j=0;j<2;++j)
			tree[tr[i]].push_back( new_triangle() );
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
		neighbor[son[0]][2] = q_ng[(ed[i]+2)%3];
	
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

	for(int i=0;i<2;++i)
	{
		vector<int>& son = tree[tr[i]];

		legalize_edge(son[0], 2);
		legalize_edge(son[1], 1);
	}
}

/** Insert a point inside a trinagle
* 
*/
void delaunay::add_point_inside(int tr_id, int pt_id)
{
	assert(tree[tr_id].empty()); // should be empty
	
	// creates 3 new triangles
	for(int i=0;i<3;++i)
		tree[tr_id].push_back( new_triangle() );
		// dosnt work with son !!!!! why ???

	vector<int>& son = tree[tr_id];
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
	for(int i=0;i<3;++i)
		legalize_edge(tree[tr_id][i], 1);
		//legalize_edge(son[i], 1); DOESN'T WORK !!!!
}
/**
*@param pt Point to be inserted to the delaunay triangulatoin
*/
void delaunay::add_point(point p)
{

	// id of the triangle containing the point
	int tr_id = search_triangle( p );
	
	int pt_id = pt.size();
	pt.push_back(p); // adds to the set of points,
			  // note that it assumes that points are unique

	int ed_id; // edge of this triangle containing the point 
	if( (ed_id = border_triangle(tr_id, p)) != -1 )
	{ // border
		add_point_border(tr_id, pt_id, ed_id);	
	}
	else
	{ // inside,
		add_point_inside(tr_id, pt_id);
	}
			
}

void delaunay::plot_points() const
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

/**
* Checks if the delaunay triangulation satisfies the delaunay condition
*O(n^2) bruteforce check
*/
bool delaunay::check() const
{
	vector< vector<int> > adj( pt.size(), vector<int>() ); // graph of points

	for(int i=0;i<triangle.size();++i)
		if(tree[i].empty()) // leaaf, actual triangle
		{
			const vector<int>& tr = triangle[i];

			for(int i=0;i<pt.size();++i)
				if(i!=tr[0] && i!=tr[1] && i!=tr[2])
					if( inside_circumcircle(pt[i], pt[tr[0]], pt[tr[1]], pt[tr[2]]) )
						return 0;
		}
	return 1;

}

void delaunay::plot_triangulation() const
{
	vector< vector<int> > adj( pt.size(), vector<int>() ); // graph of points

	for(int i=0;i<triangle.size();++i)
		if(tree[i].empty()) // leaf
		{
			for(int j=1;j<4;++j)
			{
				int a = triangle[i][j-1], b = triangle[i][j%3];
				adj[ a ].push_back( b );
			}
		}
	
	// only for the exterior triangle, adds the double link
	for(int j=1;j<4;++j)
	{
		int a = triangle[0][j-1], b = triangle[0][j%3];
		adj[ b ].push_back( a );
	}

	FILE* file = fopen("data_p","w");

	double minx,maxx,miny,maxy;

	minx = miny = 1./0.;
	maxx = maxy = -1./0.;
	//write points
	for(int i=0;i<pt.size();++i) // not artificial points
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
	fprintf(pf, "set style arrow 1 nohead\n");

	// print the graph
	vector<int> vis(pt.size(),0);

	queue<int> q;
	q.push(0);
	while(!q.empty())
	{
		int v = q.front();q.pop();
		for(int i=0;i<adj[v].size();++i)
		{
			int u=adj[v][i];
			if( vis[u] != 2) 
			{ // show edge
				fprintf(pf, "set arrow from %lf,%lf to %lf,%lf as 1\n", pt[v].x, pt[v].y, pt[u].x, pt[u].y);
			}
			//conquer
			if( vis[u] == 0 )
			{
				vis[u] = 1;
				q.push(u);
			}

		}
	}
	
	fprintf(pf, "plot \"data_p\" using 1:2 with dots title \"Points\"\n");
	fclose(pf);
	
	int result = system("gnuplot -persist graph.conf");

}

int delaunay::size() const
{
	return  triangle.size();
}
