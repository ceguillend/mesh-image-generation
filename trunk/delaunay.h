/**
 * @file Contains all the implementation of the delaunay data structure.
 */
#ifndef DELAUNAY_H_
#define DELAUNAY_H_

#include "basic_geo.h"
#include "opencv2/imgproc/imgproc.hpp"
#include <boost/functional/hash.hpp>
#include <list>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace mesh_generation {

typedef std::unordered_set< std::pair<int, int>,
        boost::hash< std::pair<int, int> > > PairSet;

typedef std::unordered_map< std::pair<int, int>, int,
        boost::hash< std::pair<int, int> > > PairMap;

/**
 * Defines the data type used to store the information of a single triangle in
 * the delaunay triangulation data structure.
 */
struct DelaunayTriangle {
  DelaunayTriangle();
  DelaunayTriangle(std::vector<int> point_ids,
                   std::vector<int> triangle_neighbor_ids);
  /**
   * Holds a triplet of point ids describing the triangle, sorted in clockwise
   * order.
   */
  std::vector<int> point_ids;

  /**
   * Holds a triplet of triangle ids describing the neighboring triangles.
   */
  std::vector<int> triangle_neighbor_ids;

  /**
   * Holds the list of triangle ids of its children (generated by splitting the
   * triangle.
   */
  std::vector<int> triangle_child_ids;
};

/**
 * Implements the delaunay triangulation algorithm described at:
 *
 * Computational Geometry: Algorithms and Applications Third Edition (March
 * 2008) Chapter 9
 *
 * As this is a offline algorithm, this implementation assumes that the points
 * are given between the rectangle [0:W,0:H]. For this particular case the
 * maximum lexicografical coordinate is always going to be p0 = (W,H), and is
 * allays assumed to be part of the set of points P (which allows the algorithm
 * to behave online).
 *
 * Note that this can be easily modified to work in offline mode, by just
 * calculating P0 on the whole set of input points.
 */
class DelaunayMesh {
public :
  /**
   * Does not initialize the delaunay mesh at all.
   */
  DelaunayMesh();

  /**
   * Creates a delaunay mesh of an image with the given height and width. it
   * tries to build the mesh trying to satisfy the segment constrains defined
   * by the simplified_cuves, and it will refine each segment up to a
   * num_split_operations times.
   */
  DelaunayMesh(int height, int width,
               const std::list< vector<Point> >& simplified_curves,
               int num_split_operations);

  /**
   * Builds the delaunay triangulation by inserting the points that appear
   * in the segment_constrains, if one of the segments doesn't show up after
   * inserting the points an assertion error is thrown.
   */
  DelaunayMesh(int height, int width,
               const std::vector< std::pair<Point,
                                            Point> >& segment_constrains);

  /**
   * Checks if the current triangulation satisfies the delaunay properties.
   * Note that this process is very expensive since it takes O(N^2) operations,
   * where N is the number of points.
   */
  bool IsValidDelaunay() const;

  /**
   * Plots the set of points already inserted. The function writes into a file
   * the commands required to plot the points using gnuplot. The files used are:
   * points_set, containing the description of the points and graph.conf,
   * containing the plotting commands.
   */
  void PlotPoints() const;

  /**
   * Plots the set of points already inserted. The function writes into a file
   * the commands required to plot the points using gnuplot. The files used are:
   * points_set, containing the description of the points and graph.conf,
   * containing the plotting commands.
   *
   * @param wrong circles Defines whether the plot should check the triangles
   *        and display the ones that don't satisfy the delaunay property, this
   *        procedure takes O(N^2) operations, where N is the number of points.
   * @param rand_circles Plots random circumcircles, plotting around 4% of the
   *        total number of circumcircles.
   * @param kill_other_plots defines if other plotting process are going to be
   *        killed or not.
   */
  void PlotTriangulation(bool wrong_circles, bool rand_circles,
                         bool display_result, bool kill_other_plots) const;

  /**
   * @return The total number of triangles in the triangle tree.
   */
  int TriangleTreeSize() const;

  /**
   * @return The number of unique points inserted in the triangulation.
   */
  int PointSetSize() const;

  /**
   * @return The number of triangles visited on every insertion operation.
   */
  double AverageInsertionCost() const;

  /**
   * @return The reference to the safe region used to build the mesh.
   */
  const cv::Mat& GetSafeRegion() const;

  /**
   * @return The set of points that belongs to a constrained delaunay edge.
   */
  void GetConstrainedPoints(vector<Point>* constrained_points) const;

  /**
   * @return A random set of safe points, with size equal to the fraction of the
   * fraction given in the input of the total number of safe points.
   */
  void GetRandSafePoints(double fraction, vector<Point>* safe_points) const;

  /**
   * @return The set of non constrained points, after applying the modified
   * Lloyd iteration, which preserves the constrained segments.
   */
  void GetUnconstrainedAdjustedPoints(vector<Point>* adjusted_points) const;

  /**
   * Inserts a safe set of points to the current delaunay triangulation. Note
   * that every point in the set should be lexicografically less than the
   * coordinate (width, height) used in the constructor and the points should
   * not break the existing constrains, otherwise the function will rise an
   * assertion error. It is expected to have  a unique set of points.
   */
  void InsertSafePoints(const vector<Point>& safe_points);

  /**
   * @return the set of constrained edges that appear in the delaunay
   * triangulation.
   */
  void GetMinimalConstrainSet(
            vector< std::pair<Point, Point> >* minimal_constrain_set) const;

 private:
  /**
   * Inserts the given point to the delaunay triangulation.
   *
   * @param point_id The point id of the point to be inserted.
   */
  void InsertPoint(int point_id);

  /**
   * Checks whether the given point is inside or in the edge of the triangle
   * defined by _triangles[triangle_id].
   * A point is considered inside a triangle if it is to the right of all the
   * triangle's edges when visiting them in clockwise order.
   */
  bool IsInsideTriangle(int triangle_id, const Point& point) const;

  /**
   * Checks whether the given points belongs to an edge of the triangle with id
   * triangle_id.
   * The parameter edge_start_pos returns the position of the start point in
   * the array of point_ids of the edge to which point belongs.
   */
  bool BelongsTriangleBorder(const Point& point, int triangle_id,
                             int* edge_start_pos) const;
  /**
   * Finds the triangle which encloses the given point.
   * @return The triangle id (a leaf in the triangle DAG). If it belongs to an
   * edge, returns any of the neighboring triangles.
   */
  int FindEnclosingTriangle(const Point& point);// const;

  /**
   * Creates a new triangle, without data.
   * @return The triangle_id of the new created triangle.
   */
  int StoreTriangle(const DelaunayTriangle& triangle);

  /**
   * Stores the given point.
   * @return The point id.
   */
  int PointId(const Point& point);

  /**
   * @return the triangle id of the new allocated triangle.
   */
  int AllocateNewTriangle();

  /**
   * @return the point id of the new allocated point.
   */
  int AllocateNewPoint();

  /**
   * @return the reference of the triangle stored with the given id.
   */
  DelaunayTriangle& GetTriangleRef(int triangle_id);

  /**
   * @return the value of the triangle stored with the given id.
   */
  const DelaunayTriangle& GetTriangleVal(int triangle_id) const;

  /**
   * @return the point stored with the given id.
   */
  Point GetPoint(int point_id) const;

  /**
   * @return whether the current point id belongs to a infinite point.
   */
  bool IsInfinitePoint(int point_id) const;

  /**
   * Checks if the circumcircle of the given triangle does not contain a point
   * inside its region, it performs a linear check in the set of points.
   */
  bool IsValidTriangle(const DelaunayTriangle& triangle) const;

  /**
   * Updates the neighbor of the given triangle. This is called when a neighbor
   * triangle is replaced by a new inserted one.
   * @param triangle_id Triangle that is going to be updated.
   * @param old_neighbor The triangle id of the old neighbor.
   * @param new_neighbor The triangle id of the new neighbor.
   */
  void UpdateNeighbor(int triangle_id, int old_neighbor, int new_neighbor);

  /**
   * Checks if the given triangle satisfies the Delaunay property at the given
   * edge. For the edge to be legal, the circumcircle of its adjacent triangles
   * should be empty. In case it is not legal, flips the edge and legalizes it.
   *
   * @param triangle_id Refers to the triangle to be checked.
   * @param side_index The triangle's side to be checked.
   */
  void LegalizeSide(int triangle_id, int side_index);

  /**
   * Adds a point to the border of a triangle defined by the triangle_id.
   *
   * @param triangle_id The id of the triangle to be modified.
   * @param point_id The id of the point to be inserted.
   * @param side_index The side in the triangle to be modified.
   */
  void InsertBorderPoint(int triangle_id, int point_id, int side_index);

  /**
   * Inserts a point inside an already existing triangle.
   *
   * @param triangle_id The id of the triangle in which the point lies.
   * @param point_id The id of the point to be inserted.
   */

  void InsertInnerPoint(int triangle_id, int point_id);
  /**
   * Finds the index in which the given neighbor_id occurs in the list of
   * neighbors of the given triangle_id
   * @param triangle_id Triangle id where the neighbor is searched.
   * @param neighbor_id Neighbor's triangle id to be found.
   */
  int NeighborSide(int triangle_id, int neighbor_id) const;

  /**
   * Writes the set of points to a file with the given name, the points are
   * write each in a separate line, with the x and y coordinated separated with
   * a tab (gnuplot format). Return the bounding box of the saved points.
   */
  void SavePointsToFile(const std::string& file_name, double* min_x,
                        double* max_x, double* min_y, double* max_y) const;
  /**
   * Plots the delaunay triangulation and the safe circles described by the
   * constrained edges.
   */
  void GetSafeRegion(bool gnu_plot, int height, int width,
      cv::Mat* safe_region) const;

  /**
   * Obtains the current set of satisfied constrains.
   */
  void GetSatisfiedConstrains(PairSet* satisfied_constrains) const;

  /**
   * @return The safe circle of segment, described by the side_index of
   * triangle with the triangle_id.
   */
   void GetSafeCircle(int triangle_id, int side_ind, Point* center,
                      double* radius) const;

  /**
   * Set of points associated with their point ids.
   */
  std::map<Point, int> _point_to_id;

  /**
   * Set of inserted points.
   */
  std::vector<Point> _points;

  /**
   * Set of border edges, pair(point_id, point_id), where the value
   * indicates the number of times the segment can be split.
   */
  PairMap _constrains;

  /**
   * Holds the set of segments with unsatisfied constrains.
   */
  PairSet  _unsatisfied_constrains;

  /**
   * Holds whether the point ids have been inserted or not.
   */
  std::unordered_set<int> _inserted_points;

  /**
   * Each position contains.
   */
  std::vector<DelaunayTriangle> _triangles;

  /**
   * Safe region.
   */
  cv::Mat _safe_region;

  /**
   * The set of movable points (safe points).
   */
  std::unordered_set<int> movable_point_ids;

  /**
   * Number of triangles visited through all the queries.
   */
  int _total_query_operations;

};
}  // namespace mesh_generation
#endif
