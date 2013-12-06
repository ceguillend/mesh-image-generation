#include"header.h"
#include"curve.h"
#include"thinning.h"
#include"basic_geo.h"
#include"delaunay.h"

using namespace std;
using namespace cv;
// source image 

//display window
char window_name[] = "Algorithm";

// output
Mat src;
Mat src_gray; // gray image
Mat gray_blur; // gray blurred image
Mat img_canny; // canny edge
Mat img_edge; // thinned edge
Mat img_group; // pixel components
Mat img_curve_px; // curve described by pixels
Mat img_curve_pt; // simplified curve 

// canny parameters
int lowThreshold = 40; // variable, 40 good plot
int const max_lowThreshold = 100;
int const ratio = 2;
int const kernel_size = 3;

// curve_parameters
int max_d = 5; // good lines
int max_len  ; // max lentgh size
int min_cv_len = 0; // show all

// Algorithm
int show_phase = 7;
const int n_phases = 7;

void write_process()
{
	imwrite("0_original.png", src);
	imwrite("1_gray.png", src_gray);
	imwrite("2_gray_blur.png", gray_blur);
	imwrite("3_canny.png", img_canny);
	imwrite("4_thinned.png", img_edge);
	imwrite("5_pixel_group.png", img_group);
	imwrite("6_pixel_curves.png", img_curve_px);
	imwrite("7_pixel_curves_simple.png", img_curve_pt);
}

void show(int, void*)
{

	int q=0;
	if(q++ == show_phase ) 
		imshow(window_name, src);		
	if(q++ == show_phase) 
		imshow(window_name, src_gray);
	if(q++ == show_phase) 
		imshow(window_name, gray_blur);

	if(q++ == show_phase) 
		imshow(window_name, img_canny);
	if(q++ == show_phase) 
		imshow(window_name, img_edge);

	if(q++ == show_phase) 
		imshow(window_name, img_group);
	if(q++ == show_phase) 
		imshow(window_name, img_curve_px);
	if(q++ == show_phase) 
		imshow(window_name, img_curve_pt);
}

/**
 * @function CannyThreshold
 * @brief Trackbar callback - Canny thresholds input with a ratio 1:3
 */
void method(int, void*)
{
	/* Find edge transformation */	

	//Mat img_canny, img_edge;
	// Reduce noise with a kernel 3x3
	blur( src_gray, gray_blur, Size(3,3) );
	// Canny detector, result img_edge
	Canny( gray_blur, img_canny, lowThreshold, lowThreshold*ratio, kernel_size );
	// thinning
	thinning_procedure(img_canny, img_edge);

	// adds the edge to the borders of the thinned image
	int nr = img_edge.rows;
	int nc = img_edge.cols;
	
	for(int i=0;i<nr;++i)	
		img_edge.at<uchar>(i, 0) = img_edge.at<uchar>(i, nc-1) = 255;

	for(int i=0;i<nc;++i)	
		img_edge.at<uchar>(0, i) = img_edge.at<uchar>(nr-1, i) = 255;

	/* Find curves */
		
	vector< vector<pii> > px_curve; // curves of pixels (pixel coordinates)
	vector< vector<point> > sg_curve; // line segments (euclidean coordinates)
//	Mat img_group; // pixel components
//	Mat img_curve_px; // curve described by pixels
//	Mat img_curve_pt; // simplified curve 

	// find pixel curves	
	find_pixel_curves(img_edge, min_cv_len, px_curve, img_group, img_curve_px);
	// fing segment curves
	simplify_curve(px_curve, img_edge.rows, img_edge.cols, max_d, max_len, sg_curve, img_curve_pt);
	//results

	// Just builds the delaunay triangulation with the segment points
	delaunay dt(src.cols, src.rows);

	set<point> pt ; // quarantee unique poinnts
	pt.insert( point( src.cols, src.rows ) );

	for(int i=0;i<sg_curve.size();++i)
		for(int j=0;j<sg_curve[i].size();++j)
		{
			if( pt.count( sg_curve[i][j] ) == 0 )
			{
				dt.add_point( sg_curve[i][j] );
				pt.insert( sg_curve[i][j] );
			}
		}
		
	//verifies dt condition
//	assert( dt.check() );
	// plot delaunay
	dt.plot_delaunay(1, 0);		
	// delaunay info
	printf("N: %d  NT: %d  NLO: %.2lf\n", dt.npoints(), dt.size(), dt.average_location_operations());
	// show the image
	show(0,0);
	// write process
//	write_process();
}



/** @function main */
int main( int argc, char** argv )
{
	
	if( argc != 2)
	{
		printf("./method [file_name]\n");
		return -1; 
	}
	/// Load an image
	src = imread( argv[1] );

	if( !src.data )
	{ 
		printf("./method [file_name]\n");
		return -1; 
	}

	/// Convert the image to grayscale
	cvtColor( src, src_gray, CV_BGR2GRAY );

	/// Create a window
	namedWindow( window_name, CV_WINDOW_AUTOSIZE );

	int ssz = sqrt(sqr(src.rows)+sqr(src.cols));

	/// Create a Trackbar for user to enter threshold
	createTrackbar( "Phase", window_name, &show_phase, n_phases, show );
	
	/// Create a Trackbar for user to enter threshold
	createTrackbar( "3) Min Thres", window_name, &lowThreshold, max_lowThreshold, method );
	
	/// Create a Trackbar for user to enter min_cmp_sz , the minimum component size
	createTrackbar( "6) Min cv_size", window_name, &min_cv_len, ssz, method );
	
	max_len = ssz; // show all	
	/// Create a Trackbar for user to enter max_sz
	createTrackbar( "7) Max Leng", window_name, &max_len, ssz, method );
	
	/// Create a Trackbar for user to enter min_cmp_sz , the minimum component size
	createTrackbar( "7) Max Dist", window_name, &max_d, 20, method );
	
	
	/// Show the image
	method(0, 0);

	/// Wait until user exit program by pressing a key
	waitKey(0);

	return 0;
}
