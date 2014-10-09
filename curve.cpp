#include "curve.h"

using namespace std;
using namespace cv;

const int MAX_LUM = 200;
void find_pixel_curves(const Mat& img_edge, int min_curve_len,
                        vector< vector<pii> >& px_curve, Mat& img_group,
                        Mat& img_curve  ) {
	assert( img_edge.type() == CV_8UC1 );
	srand(SEED);
	// 8-neighborhood
	const int r_v[] = {-1,1,0,0,-1, 1,1,-1};
	const int c_v[] = {0,0,-1,1, 1,-1,1,-1};
	
	// dimension
	const int n = img_edge.rows;
	const int m = img_edge.cols;

	// image of groups
	img_group = Mat::zeros(n,m,CV_8UC3);

	// defines the curve to which the given pixels belongs
	vector< vector<int> > id_cv(n, vector<int>(m, -1)); 
	int l_curve = -1 ;

	vector< vector<pii> > temp_curve; 
	// finding groups and curves
	for(int i=0;i<n;++i)
		for(int j=0;j<m;++j)
			if( img_edge.at<uchar>(i, j) == 255 && id_cv[i][j]==-1 )
			{
				Vec3b group_color = rand_color();
				
				queue< pii > q;
				id_cv[i][j] = ++l_curve;
				temp_curve.push_back( vector<pii>() );
				q.push( pii(i,j) );		

        // bfs
				while(!q.empty())
				{
					int r = q.front().F, c = q.front().S;
					q.pop();
					img_group.at<Vec3b>(r, c) = group_color;
					temp_curve[id_cv[r][c]].push_back( pii(r, c) ); 

					//expand
					for(int d=0, cnt=0;d<8;++d)
					{
						int nr = r+r_v[d], nc = c+c_v[d];
						if( 0 <= nr && nr < n && 0 <= nc && nc < m ) 
						{
							if( img_edge.at<uchar>(nr, nc) == 255 && id_cv[nr][nc]==-1 )
							{
								if( cnt )
								{ // split
									id_cv[nr][nc] = ++l_curve;
									temp_curve.push_back( vector<pii>() );
									// to maintain connectivity
									//temp_curve.back().push_back( pii(r,c) );
								}
								else
								{
									id_cv[nr][nc] = id_cv[r][c];	
								}
								q.push( pii(nr,nc) );
								++cnt;
							}
						}
					}
				}
			}

	/* ploting curves */

	px_curve.clear();
	
	//image curve plottinh	
  // white background
	img_curve = Mat::zeros(n,m,CV_8UC3);
	srand(SEED);	
	for(int i=0;i<temp_curve.size();++i)
	{
		if( temp_curve[i].size() >= min_curve_len )
		{
			px_curve.push_back( temp_curve[i] );

			Vec3b curve_color = rand_color();
			for(int j=0;j<temp_curve[i].size();++j)
			{
				img_curve.at<Vec3b>(temp_curve[i][j].F, temp_curve[i][j].S) =
            curve_color;
			}
		}
	}
}


void simplify_curve(vector< vector<pii> >& px_curve, int H, int W, double max_d,
                    double max_len, vector< vector<point> >& sg_curve,
                    Mat& img_curve_pt) {
	sg_curve.assign( px_curve.size(),  vector<point>() );

	for(int i=0;i<px_curve.size();++i)
	{
		vector<point> seq(px_curve[i].size());

		for(int j=0;j<seq.size();++j)
			seq[j] = px_pt(px_curve[i][j], H);
		
		stack< pii > st;
		st.push( pii(0, seq.size()-1) );

		while(!st.empty())
		{
			int l = st.top().F, r = st.top().S;
			st.pop();
			
			bool split = cmp( seq[l].dist(seq[r]), max_len) > 0; 
			int pos = (l+r)/2;

			if( !split )
			{
				// find equation
				double a = seq[r].y-seq[l].y;
				double b = -(seq[r].x-seq[l].x);
				double c = -( a*seq[l].x + b*seq[l].y );
				double norm = sqrt( sqr(a)+sqr(b) );

				double dist=0;
				//find farthest point
				for(int j = l+1; j < r; ++j)
				{
					double qd = fabs(a*seq[j].x+b*seq[j].y+c)/norm;	
					if( cmp(dist, qd ) < 0 )
					{
						dist = qd;
						pos = j;
					}
				}
				split = cmp(dist, max_d ) > 0 ;
			}

			if( split )
			{
				st.push( pii(pos,r) );
				st.push( pii(l,pos) );//left last to visit firs
			}
			else // set this
			{
				sg_curve[i].push_back( seq[l] ); // only first to not repeat points
			}
		}
		sg_curve[i].push_back( seq.back() );
	}
	//image curve	
	img_curve_pt = Mat::zeros(H,W,CV_8UC3);
	
	srand(SEED);	
	for(int i=0;i<sg_curve.size();++i)
	{
    Vec3b rgb_color = rand_color();
		Scalar curve_color(rgb_color[0], rgb_color[1], rgb_color[2]);
		pii px, ppx;
		for(int j=0;j<sg_curve[i].size();++j)
		{
			px = pt_px( sg_curve[i][j], H );
			circle(img_curve_pt, cvPoint(px.S, px.F), 3, curve_color);
			if(j)
			{
				line(img_curve_pt, cvPoint(ppx.S,ppx.F), cvPoint(px.S,px.F),
              curve_color);
			}
			ppx = px;
		}
	}

}
