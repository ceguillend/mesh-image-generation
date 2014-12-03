/** @file
 *  Executes the mesh image generation algorithm.
 */

#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "curve.h"
#include "thinning.h"
#include "basic_geo.h"
#include "delaunay.h"
#include "time.h"
#include <set>
#include <list>
#include <stdio.h>

using namespace std;
using cv::Mat;
using cv::blur;
using cv::Canny;
using cv::Size;
using cv::waitKey;
using cv::imread;
using cv::createTrackbar;
using std::list;
using std::set;
using std::vector;
/**
 * OpenCv Display window."
 */
const char kWindowResult[] = "Image Result";
/**
 * Canny Parameters
 */
int low_threshold = 40;
const int kMaxLowThreshold = 100;
const int kRatio = 2;
const int kKernelSize = 3;

/**
 * Curve parameters.
 */
int component_size_threshold = 0;
int curve_length_threshold = 0;
int max_point_distance = 1;
int max_segment_length;
const int kMaxPointDistance = 30;

/**
 * Algorithm parameters.
 */
int display_phase = 0;
const int kNumPhases = 6; // 0 .. kNumPhases
// Involves a O(N^2) check.
const bool kPlotWrongCircumcircle = false;
// Samples 4 percent of the triangles.
const bool kPlotRandCircumcircle = false;
// Defines whether it generates the plot or not.
const bool kDisplayTriangulation = true;
// Defines whether a gnu plot is generated with with the safe circles.
const bool kGnuDisplaySafeCircles = true;
// Number of fixing splits.
const int kNumSplitOperations = 5;

namespace mesh_generation {

/**
 * Source image.
 */
Mat image_source;

/**
 * Gray scale image.
 */
Mat image_gray;

/**
 * Gray scale blurred image.
 */
Mat image_gray_blur;

/**
 * Canny edge detected.
 */
Mat image_canny;

/**
 * Thinned 1-edge.
 */
Mat image_thinned_edge;

/**
 * Pixel based curves found.
 */
Mat image_pixel_curve;

/**
 * Simplified point curve.
 */
Mat image_point_curve;

/**
 * Curves used to generate the delaunay triangulation.
 */
list< vector<Point> > simplified_curves;


void write_process()
{
	imwrite("0_original.png", image_source);
	imwrite("1_gray.png", image_gray);
	imwrite("2_gray_blur.png", image_gray_blur);
	imwrite("3_canny.png", image_canny);
	imwrite("4_thinned.png", image_thinned_edge);
	imwrite("5_pixel_curves.png", image_pixel_curve);
	imwrite("6_pixel_curves_simple.png", image_point_curve);
}

void DisplayImageResult(int, void*) {
	int q=0;
	if (q++ == display_phase) {
		imshow(kWindowResult, image_source);
  }
	if (q++ == display_phase) {
		imshow(kWindowResult, image_gray);
  }
	if (q++ == display_phase) {
		imshow(kWindowResult, image_gray_blur);
  }
	if (q++ == display_phase) {
    imshow(kWindowResult, image_canny);
  }
	if (q++ == display_phase) {
		imshow(kWindowResult, image_thinned_edge);
  }
	if (q++ == display_phase) {
		imshow(kWindowResult, image_pixel_curve);
  }
	if (q++ == display_phase) {
		imshow(kWindowResult, image_point_curve);
  }
}

/**
 * Displays the stats of the current execution.
 * @param method_time Time elapsed in seconds during the method execution.
 * @param dt the delaunay triangulation generated.
 * @param plot_time Time elapsed in seconds during the plot execution.
 */
void PrintTime(const string& method_name, int method_time) {
  printf("====================================================\n");
  printf("%s:\n", method_name.c_str());
  printf("  -Time elapsed: %3d h %2d m %2d s\n\n", method_time/60/60,
          (method_time/60)%60, method_time%60);
}

/**
 * Runs the mesh image generation algorithm.
 */
void CurvesGeneration(int, void*) {
  clock_t begin_time = clock();
	int num_rows = image_source.rows;
	int num_cols = image_source.cols;
	// Converts the image to grayscale.
	cvtColor(image_source, image_gray, CV_BGR2GRAY );
	// Reduce noise with a kernel 3x3.
	blur(image_gray, image_gray_blur, Size(3, 3));
  // Finds the edges.
	Canny(image_gray_blur, image_canny, low_threshold, low_threshold * kRatio,
          kKernelSize);
  // Thins the edges.
	ThinningProcedure(image_canny, &image_thinned_edge);
	// Adds edge to the borders of the thinned image to fully cover the image
  // domain.
  for (int i = 0; i < num_rows; ++i) {
		image_thinned_edge.at<uchar>(i, 0) =
        image_thinned_edge.at<uchar>(i, num_cols - 1) = 255;
  }

  for (int  i = 0; i < num_cols; ++i) {
		image_thinned_edge.at<uchar>(0, i) =
        image_thinned_edge.at<uchar>(num_rows - 1, i) = 255;
  }

	list< vector<Point> > point_curves;
	FindPixelCurves(image_thinned_edge, component_size_threshold,
                    curve_length_threshold, &point_curves, &image_pixel_curve);

  // Sets the global variable of simplified_curves, to be used in the mesh
  // construction.
	SimplifyPointCurves(point_curves, num_rows, num_cols, max_point_distance,
                        max_segment_length, &simplified_curves,
                        &image_point_curve);

	DisplayImageResult(0, nullptr);

  clock_t end_time = clock();
  PrintTime("Image Processing", (end_time-begin_time) / CLOCKS_PER_SEC);
  // write_process();
}
// Generates the delaunay triangulation.

void PrintDenaulayStats (const DelaunayMesh& mesh) {
  printf("Delaunay:\n");
	printf("  - Points inserted     : %d\n", mesh.PointSetSize());
  printf("  - Triangle tree nodes : %d\n", mesh.TriangleTreeSize());
  printf("  - Search point cost   : %.2lf\n", mesh.AverageInsertionCost());
}

void MeshGeneration(int, void*) {
  clock_t begin_time = clock();
  DelaunayMesh mesh(image_source.rows, image_source.cols, simplified_curves,
                    kNumSplitOperations);
	// Verifies mesh correctness O(n^2).
  clock_t end_time = clock();
  assert(mesh.IsValidDelaunay());
  mesh.PlotTriangulation(false, false, true);

  PrintTime("Delaunay Triangulation", (end_time-begin_time) / CLOCKS_PER_SEC);
  PrintDenaulayStats(mesh);
  imshow(kWindowResult, mesh.GetSafeRegion());
}
}  // namespace mesh_generation

int main( int argc, char** argv ) {

	if ( argc != 2) {
		printf("./method [file_name]\n");
		return -1;
	}
	// Loads the image.
	mesh_generation::image_source = imread(argv[1]);

	if (!mesh_generation::image_source.data) {
		printf("./method [file_name]\n");
		return -1;
	}

	// Creates a cv window.
  cv::namedWindow(kWindowResult, CV_WINDOW_AUTOSIZE);

	// Creates a trackbar for user to enter the phase to be displayed.
	createTrackbar("Phase", "", &display_phase, kNumPhases,
                    mesh_generation::DisplayImageResult);

	// Creates a trackbar for user to enter threshold
	createTrackbar("Canny Threshold", "", &low_threshold, kMaxLowThreshold,
                    mesh_generation::CurvesGeneration);

  max_segment_length =
      sqrt(mesh_generation::Sqr(mesh_generation::image_source.rows) +
            mesh_generation::Sqr(mesh_generation::image_source.cols));

  // Creates a trackbar allow the definition of the minimum allowed pixel
  // component size.
	createTrackbar("Component Size Threshold", "",
                    &component_size_threshold,
                    max_segment_length, mesh_generation::CurvesGeneration);

  // Creates a trackbar to define the minimum curve length allowed.
	createTrackbar("Min Curve Length", "",
                    &curve_length_threshold, max_segment_length,
                    mesh_generation::CurvesGeneration);


	// Creates a trackbar for user to enter the maximum allowed segment length.
	createTrackbar("Max Segment Length", "", &max_segment_length,
                    max_segment_length, mesh_generation::CurvesGeneration);

	// Creates a trackbar for user to enter the maximum allowed distance to the
  // simplified curve.
	createTrackbar("Max Point Distance", "", &max_point_distance,
                     kMaxPointDistance, mesh_generation::CurvesGeneration);

  // Runs the algorithm, withe the paremets
  cv::createButton("Mesh Processing", mesh_generation::MeshGeneration);

	// Show the image
	mesh_generation::CurvesGeneration(0, nullptr);
	// Wait until user exit program by pressing the ESC key
  while(cv::waitKey(1) != 27);

	return 0;
}
