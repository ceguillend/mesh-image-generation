#include "delaunay.h"
#include <assert.h>
#include <stdio.h>
#include <algorithm>

using std::vector;
using std::swap;
using std::max;
using std::min;
using std::string;
using std::pair;
using std::make_pair;
using std::list;
using cv::Mat;

const int kTriangleSides = 3;
const double kInfinite = 1./0.;
const int kP2Id = 0;
const int kP1Id = 1;
const int kNullId = -1;
const int kTriangleRoot = 0;
const bool kPlotSafeRegion = false;
const int kSafeRegion = 255;

namespace {
pair<int, int> minMaxPair(int a, int b) {
  return make_pair(min(a, b), max(a,b));
}
}  // namespace

namespace mesh_generation {

DelaunayTriangle::DelaunayTriangle():
  point_ids(kTriangleSides, kNullId),
  triangle_neighbor_ids(kTriangleSides, kNullId) {}

DelaunayTriangle::DelaunayTriangle(std::vector<int> point_ids,
                                   std::vector<int> triangle_neighbor_ids):
  point_ids(point_ids),
  triangle_neighbor_ids(triangle_neighbor_ids) {
  assert(point_ids.size() == 3);
  assert(triangle_neighbor_ids.size() == 3);
}

DelaunayMesh::DelaunayMesh(int height, int width,
                           const list< vector<Point> >& simplified_curves,
                           int num_split_operations) {
  assert(num_split_operations >= 0);
  _total_query_operations = 0;
  Point lexicographically_bigger(width, height);
  // Creates the "outside-triangle" which is big enough to contain any triangle
  // inserted, the coordinates of the imaginary points are not actually used.
  int p2_id = PointId(Point(kInfinite, kInfinite));
  int p1_id = PointId(Point(-kInfinite, -kInfinite));
  int p0_id = PointId(lexicographically_bigger);

  int outer_triangle_id = AllocateNewTriangle();
  GetTriangleRef(outer_triangle_id).point_ids = {p2_id, p0_id, p1_id};

  // Stores the points and constrains.
  for (const vector<Point>& curve: simplified_curves) {
    assert(curve.size() > 1);
    for (int i = 0; i < curve.size(); ++i) {
      assert(Cmp(0, curve[i].x) <= 0 && Cmp(curve[i].x, width-1) <= 0);
      assert(Cmp(0, curve[i].y) <= 0 && Cmp(curve[i].y, height-1) <= 0);
      if (i) {
        assert(curve[i] != curve[i - 1]);
        int point1 = PointId(curve[i]);
        int point2 = PointId(curve[i - 1]);
        assert(point1 != point2);
        _constrains[minMaxPair(point1, point2)] = num_split_operations;
        _unsatisfied_constrains.insert(minMaxPair(point1, point2));
      }
    }
  }
  vector<int> point_ids;
  // Shuffles the points.
  for (int point_id = 3; point_id < _points.size(); ++point_id) {
    point_ids.push_back(point_id);
  }
  for (int i = 0; i < point_ids.size(); ++i) {
    swap(point_ids[i], point_ids[rand()%(i+1)]);
  }
  // Inserts the points in the triangulation.
  for (int i = 0; i < point_ids.size(); ++i) {
    InsertPoint(point_ids[i]);
  }
  printf("\nInitial unsatisfied constrains: %d\n",
          (int)_unsatisfied_constrains.size());
  int num_inserted_points = 0;
  // Tries to solve the unsatisfied constrains.
  while(!_unsatisfied_constrains.empty()) {
    auto it = _unsatisfied_constrains.begin();
    assert(!IsInfinitePoint(it->first) && !IsInfinitePoint(it->second));
    int new_point = PointId((GetPoint(it->first) + GetPoint(it->second)) / 2);
    _constrains[minMaxPair(new_point, it->first)] = _constrains[*it] - 1;
    _constrains[minMaxPair(new_point, it->second)] = _constrains[*it] - 1;
    _constrains[*it] = 0;
    _unsatisfied_constrains.erase(it);
    //  Inserts the point to the delaunay triangulation.
    InsertPoint(new_point);
    ++num_inserted_points;
  }
  printf("Number of additional points inserted: %d\n", num_inserted_points);
  PairSet initial_constrains;
  PairSet final_constrains;

  GetSatisfiedConstrains(&initial_constrains);
  // Obtains the safe region
  GetSafeRegion(kPlotSafeRegion, height, width, &_safe_region);
  vector<Point> inner_region_points;
  for (int row = 0; row < height; ++row) {
    for (int col = 0; col < width; ++col) {
      if (_safe_region.at<uchar>(row, col) != kSafeRegion) {
        inner_region_points.push_back(Pixel(row, col).ToPoint(height));
      }
    }
  }
  int sample_size = min(PointSetSize() / 3,
                        (int)inner_region_points.size()/100);
  for (int i = 0; i < inner_region_points.size(); ++i) {
    int p = rand() % (i + 1);
    if (p < sample_size) {
      swap(inner_region_points[p], inner_region_points[i]);
    }
  }
  for (int i = 0; i < sample_size; ++i) {
    InsertPoint(PointId(inner_region_points[i]));
  }
  printf("Number of inner points inserted: %d\n", num_inserted_points);
  GetSatisfiedConstrains(&final_constrains);

  // Checks the constrains preservation.
  for (const pair<int, int>& constrain: initial_constrains) {
    assert(final_constrains.count(constrain) > 0);
  }
}

int DelaunayMesh::PointId(const Point& point) {
  if (_point_to_id.count(point)) {
    return _point_to_id[point];
  }

  int point_id = _points.size();
  _points.push_back(point);
  return _point_to_id[point] = point_id;
}

Point DelaunayMesh::GetPoint(int point_id) const {
  return _points[point_id];
}

int DelaunayMesh::AllocateNewTriangle() {
  _triangles.push_back(DelaunayTriangle());
  return _triangles.size() - 1;
}

DelaunayTriangle& DelaunayMesh::GetTriangleRef(int triangle_id) {
  return _triangles[triangle_id];
}

const DelaunayTriangle& DelaunayMesh::GetTriangleVal(int triangle_id) const {
  return _triangles[triangle_id];
}

bool DelaunayMesh::IsInfinitePoint(int point_id) const {
  // Ids 0 (p2) and 1(p1) correspond to the abstraction of infinite points.
  return point_id == kP2Id || point_id == kP1Id;
}


bool DelaunayMesh::IsInsideTriangle(int triangle_id, const Point& point) const {
  for (int pos = 0; pos < kTriangleSides; ++pos) {
    const int beg_point_id = GetTriangleVal(triangle_id).point_ids[pos];
    const int end_point_id = GetTriangleVal(triangle_id)
                                  .point_ids[(pos + 1)% 3];

    // The edges goes from beg -> end, in clockwise order.
    const Point beg = GetPoint(beg_point_id);
    const Point end = GetPoint(end_point_id);

    bool beg_inf = IsInfinitePoint(beg_point_id);
    bool end_inf = IsInfinitePoint(end_point_id);

    if (!beg_inf && !end_inf) {
      if (Cmp((end - beg) ^ (point - beg), 0) > 0) {
        return false;
      }
    } else if (beg_inf && !end_inf) {
      if (beg_point_id == kP2Id && point > end) {
        return false;
      } else if (beg_point_id == kP1Id && point < end) {
        return false;
      }
    } else if (!beg_inf && end_inf) {
      if (end_point_id == kP2Id && point < beg) {
        return false;
      } else if (end_point_id == kP1Id && point > beg) {
        return false;
      }
    }  // else: beg_inf && end_inf, the point is considered always inside.
  }
  return true;
}


bool DelaunayMesh::BelongsTriangleBorder(const Point& point, int triangle_id,
                                         int* edge_start_pos) const {
  for (int pos = 0; pos < kTriangleSides; ++pos) {
    int beg_point_id = GetTriangleVal(triangle_id).point_ids[pos];
    int end_point_id = GetTriangleVal(triangle_id).point_ids[(pos + 1)% 3];

    // if one of the points are the infinite points, then the point does not
    // belong to it by definition). There is no point co-linear with p1 and p2
    if (!IsInfinitePoint(beg_point_id) && !IsInfinitePoint(end_point_id)) {
    const Point& beg = GetPoint(beg_point_id);
      const Point& end = GetPoint(end_point_id);
      if (Cmp((end - beg)^(point - beg), 0) == 0) {
        *edge_start_pos = pos;
        return true;
      }
    }
  }
  return false;
}

int DelaunayMesh::FindEnclosingTriangle(const Point& point) {
  int triangle_id = kTriangleRoot; // root

  while (!GetTriangleVal(triangle_id).triangle_child_ids.empty()) {
    int next_triangle_id = kNullId;
    for (int child_triangle_id :
              GetTriangleVal(triangle_id).triangle_child_ids) {
      ++_total_query_operations;
      if (IsInsideTriangle(child_triangle_id, point)) {
        next_triangle_id = child_triangle_id;
        break;
      }
    }
    assert (next_triangle_id != kNullId);
    triangle_id = next_triangle_id;
  }
  return triangle_id;
}


void DelaunayMesh::UpdateNeighbor(int triangle_id, int old_neighbor,
                                  int new_neighbor) {
  if (triangle_id == kNullId) {
    return;
  }
  for (int& neighbor_id : GetTriangleRef(triangle_id).triangle_neighbor_ids) {
    if (neighbor_id == old_neighbor) {
      neighbor_id = new_neighbor;
      return;
    }
  }
  assert(false);
}


int DelaunayMesh::NeighborSide(int triangle_id, int neighbor_id) const {
  for (int side_index = 0; side_index < kTriangleSides; ++side_index) {
    if (neighbor_id ==
            GetTriangleVal(triangle_id).triangle_neighbor_ids[side_index]) {
      return side_index;
    }
  }
  // It is only called when the answer exists
  assert(false);
}


void DelaunayMesh::LegalizeSide(int triangle_id, int side_index) {
  int neighbor_id = GetTriangleVal(triangle_id)
                        .triangle_neighbor_ids[side_index];
  if (neighbor_id == kNullId) {
    // Checks the case of the outside triangle, always valid.
    return ;
  }
  int neighbor_side_index = NeighborSide(neighbor_id, triangle_id);
  // Related points.
  // At most one point in the side can be infinite.
  int I = GetTriangleVal(triangle_id).point_ids[side_index];
  int J = GetTriangleVal(triangle_id)
              .point_ids[(side_index + 1) % kTriangleSides];
  int K = GetTriangleVal(triangle_id)
              .point_ids[(side_index + 2) % kTriangleSides];
  int L = GetTriangleVal(neighbor_id)
              .point_ids[(neighbor_side_index + 2) % kTriangleSides];
  // K is check point (initially the newly inserted point), never infinite.
  assert(!IsInfinitePoint(K));
  bool ilegal_edge = false;
  if (!IsInfinitePoint(I) &&
      !IsInfinitePoint(J) &&
      !IsInfinitePoint(L)) {
    // All are normal points.
    ilegal_edge = IsInsideCircumcircle(GetPoint(L),
                                       GetPoint(I), GetPoint(J), GetPoint(K));
  } else {
    if(!IsInfinitePoint(L)) {
      // One of the elements of the edge is an infinite point.
      if(GetPoint(L) > GetPoint(K) ) {
        swap(L, K);
      }
      // GetPoint(L) < GetPoint(K), lexicografical orter.
      if (I > J) {
        swap(I, J);
      }
      // I: holds the infinite vertex.
      assert(GetPoint(L) < GetPoint(J) && GetPoint(J) < GetPoint(K));

      if (I == kP1Id) { // p-1
        // Flips only if convex.
        ilegal_edge = Cmp((GetPoint(J) - GetPoint(L)) ^
                          (GetPoint(K) - GetPoint(L)) , 0) < 0 ;
      } else if (I == kP2Id) { // I == 0 // p-2
        // Flips only if convex.
        ilegal_edge = Cmp((GetPoint(J) - GetPoint(L)) ^
                          (GetPoint(K) - GetPoint(L)) , 0) > 0 ;
      }
    } else {
      ilegal_edge = false;
    }
  }

  if (ilegal_edge) {
    vector<int> split_triangles = {triangle_id, neighbor_id};
    vector<int> side_indices = {side_index, neighbor_side_index};
    vector<int> triangle_children = {AllocateNewTriangle(),
                                      AllocateNewTriangle()};
    // Adjust the neighbor values and points of the new triangles.
    for (int i = 0; i < triangle_children.size(); ++i) {
      const auto& current_triangle = GetTriangleVal(split_triangles[i]);
      const auto& oposite_triangle = GetTriangleVal(split_triangles[i ^ 1]);
      const int triangle_child = triangle_children[i];

      GetTriangleRef(triangle_child).point_ids = {
        current_triangle.point_ids[(side_indices[i] + 2) % 3],
        current_triangle.point_ids[side_indices[i]],
        oposite_triangle.point_ids[(side_indices[i ^ 1] + 2) % 3]
      };

      GetTriangleRef(triangle_child).triangle_neighbor_ids = {
        current_triangle.triangle_neighbor_ids[(side_indices[i] + 2) % 3],
        oposite_triangle.triangle_neighbor_ids[(side_indices[i ^ 1] + 1) % 3],
        triangle_children[i ^ 1]
      };

      UpdateNeighbor(GetTriangleVal(triangle_child).triangle_neighbor_ids[0],
                     split_triangles[i], triangle_child);
      UpdateNeighbor(GetTriangleVal(triangle_child).triangle_neighbor_ids[1],
                     split_triangles[i ^ 1], triangle_child);
    }
    // Adds children to the split triangles.
    for (int split_triangle: split_triangles) {
      GetTriangleRef(split_triangle).triangle_child_ids = triangle_children;
    }
    // Updates the constrains information.
    _unsatisfied_constrains.erase(minMaxPair(K, L));
    const auto it = _constrains.find(minMaxPair(I, J));
    if (it != _constrains.end() && it->second > 0) {
      if (!IsInfinitePoint(I) && !IsInfinitePoint(J)) {
        _unsatisfied_constrains.insert(minMaxPair(I, J));
      }
    }
    // Expands legalization.
    LegalizeSide(triangle_children[0], 1);
    LegalizeSide(triangle_children[1], 0);
  }
}

void DelaunayMesh::InsertBorderPoint(int triangle_id, int point_id,
                                     int side_index) {
  // Checks that the triangle is a leaf.
  assert(GetTriangleVal(triangle_id).triangle_child_ids.empty());
  // Gets the triangle neighbor.
  int neighbor_id = GetTriangleVal(triangle_id)
                        .triangle_neighbor_ids[side_index];
  // Checks that the triangle is a leaf.
  assert(GetTriangleVal(neighbor_id).triangle_child_ids.empty());
  // Neighbor's side.
  int neighbor_side_index = NeighborSide(neighbor_id, triangle_id);

  vector<int> split_triangles = {triangle_id, neighbor_id};
  vector<int> side_indices = {side_index, neighbor_side_index};

  // Creates two children per triangle.
  vector<int> triangle_children[] = {
       {AllocateNewTriangle(), AllocateNewTriangle()},
       {AllocateNewTriangle(), AllocateNewTriangle()}
  };

  // Updates the neighborhood and sets the new triangle's points.
  for (int i = 0; i < split_triangles.size(); ++i) {
    const auto& current_triangle = GetTriangleVal(split_triangles[i]);
    const int current_side = side_indices[i];
    vector<int>& current_children = triangle_children[i];
    vector<int>& oposite_children = triangle_children[i ^ 1];

    // First child.
    GetTriangleRef(current_children[0]).point_ids = {
      current_triangle.point_ids[current_side],
      point_id,
      current_triangle.point_ids[(current_side + 2) % 3]
    };

    GetTriangleRef(current_children[0]).triangle_neighbor_ids = {
      oposite_children[1],
      current_children[1],
      current_triangle.triangle_neighbor_ids[(current_side + 2) % 3]
    };

    UpdateNeighbor(GetTriangleVal(current_children[0]).triangle_neighbor_ids[2],
                   split_triangles[i], current_children[0]);

    // Second child.
    GetTriangleRef(current_children[1]).point_ids = {
      point_id,
      current_triangle.point_ids[(current_side + 1) % 3],
      current_triangle.point_ids[(current_side + 2) % 3]
    };

    GetTriangleRef(current_children[1]).triangle_neighbor_ids = {
      oposite_children[0],
      current_triangle.triangle_neighbor_ids[(current_side + 1) % 3],
      current_children[0]
    };

    UpdateNeighbor(GetTriangleVal(current_children[1]).triangle_neighbor_ids[1],
                   split_triangles[i], current_children[1]);

    // Updates constrains.
    _unsatisfied_constrains
        .erase(minMaxPair(point_id,
                          current_triangle.point_ids[current_side]));
    _unsatisfied_constrains
        .erase(minMaxPair(point_id,
                          current_triangle.point_ids[(current_side + 2) % 3]));
  }
  vector<int> edge_points = {
    GetTriangleVal(triangle_id).point_ids[side_index],
    GetTriangleVal(triangle_id).point_ids[(side_index + 1) % 3]
  };

  // Check removed edge.
  const pair<int, int> rm_edge = minMaxPair(edge_points[0], edge_points[1]);
  const auto it = _constrains.find(rm_edge);
  if (it != _constrains.end()) {
    if (it->second > 0) {
      // replaces the constrain.
      for (int i = 0; i < edge_points.size(); ++i) {
        const pair<int, int> new_edge = minMaxPair(point_id, edge_points[i]);
        _constrains[new_edge] = max(_constrains[new_edge], it->second);
      }
    }
    it->second = 0;
  }

  // Adds children.
  for (int i = 0; i < split_triangles.size(); ++i) {
    GetTriangleRef(split_triangles[i]).triangle_child_ids =
          triangle_children[i];
  }

  // Legalize edges.
  for (int i = 0; i < split_triangles.size(); ++i) {
    LegalizeSide(triangle_children[i][0], 2);
    LegalizeSide(triangle_children[i][1], 1);
  }
}

void DelaunayMesh::InsertInnerPoint(int triangle_id, int point_id) {
  // Checks that the triangle is a leaf.
  assert(GetTriangleVal(triangle_id).triangle_child_ids.empty());

  vector<int> triangle_children = {AllocateNewTriangle(), AllocateNewTriangle(),
                                   AllocateNewTriangle()};
  // Sets the appropriate values to the newly created triangles.
  for (int i = 0 ; i < triangle_children.size(); ++i) {
    const int triangle_child = triangle_children[i];

    GetTriangleRef(triangle_child).point_ids = {
      point_id,
      GetTriangleVal(triangle_id).point_ids[i],
      GetTriangleVal(triangle_id).point_ids[(i + 1) % 3]
    };

    GetTriangleRef(triangle_child).triangle_neighbor_ids = {
      triangle_children[(i + 2) % 3],
      GetTriangleVal(triangle_id).triangle_neighbor_ids[i],
      triangle_children[(i + 1) % 3]
    };
    // update the neighbor's information
    UpdateNeighbor(GetTriangleVal(triangle_id).triangle_neighbor_ids[i],
                   triangle_id, triangle_child);
  }
  // Sets the triangle's children
  GetTriangleRef(triangle_id).triangle_child_ids = triangle_children;

  // Updates constrains.
  for (int i = 0; i < kTriangleSides; ++i) {
    _unsatisfied_constrains
        .erase(minMaxPair(point_id, GetTriangleVal(triangle_id).point_ids[i]));
  }
  // Legalize edges.
  for (int triangle_child: triangle_children) {
    LegalizeSide(triangle_child, 1);
  }
}

void DelaunayMesh::InsertPoint(int point_id) {
  if (_inserted_points.count(point_id)) {
    return;
  }
  _inserted_points.insert(point_id);
  // Triangle enclosing the point to be inserted.
  int triangle_id = FindEnclosingTriangle(GetPoint(point_id));
  // Checks whether the point belongs to the inner region or an edge.
  int edge_id;
  if (BelongsTriangleBorder(GetPoint(point_id), triangle_id, &edge_id)) {
    InsertBorderPoint(triangle_id, point_id, edge_id);
  } else {
    InsertInnerPoint(triangle_id, point_id);
  }
}

void DelaunayMesh::GetSatisfiedConstrains(PairSet* satisfied_constrains) const {
  for (const DelaunayTriangle& triangle: _triangles) {
    // Plots only the leaves.
    if(triangle.triangle_child_ids.empty()) {
      // Plots connectivity among points.
      for (int i = 0; i < kTriangleSides; ++i) {
        int beg_point_id = triangle.point_ids[i];
        int end_point_id = triangle.point_ids[(i + 1) % 3];
        if (!IsInfinitePoint(beg_point_id) &&
            !IsInfinitePoint(end_point_id)) {
          const Point& beg_point = GetPoint(beg_point_id);
          const Point& end_point = GetPoint(end_point_id);
          if (beg_point == end_point) {
            fprintf(stderr, "%.6f %.6f\n%.6f %.6f\n\n", beg_point.x,
                beg_point.y, end_point.x, end_point.y);
            fprintf(stderr, "%d %d\n", beg_point_id, end_point_id);
          }
          assert (beg_point != end_point);
          const auto segment = minMaxPair(beg_point_id, end_point_id);
          // Makes sure that the edge is printed only once.
          if (beg_point_id < end_point_id && _constrains.count(segment)) {
            satisfied_constrains->insert(segment);
          }
        }
      }
    }
  }
}

void DelaunayMesh::SavePointsToFile(const string& file_name, double* min_x,
                                    double* max_x, double* min_y,
                                    double* max_y) const {
  *min_x = *min_y = 1. / 0.;
  *max_x = *max_y = -1. / 0.;
  FILE* points_file = fopen(file_name.c_str(), "w");
  // Writes the set of points.
  for (int point_id = 0; point_id < _points.size(); ++point_id) {
    if (!IsInfinitePoint(point_id)) {
      const Point& point  = _points[point_id];
      fprintf(points_file, "%lf\t%lf\n", point.x, point.y);
      *min_x = min(*min_x, point.x);
      *max_x = max(*max_x, point.x);
      *min_y = min(*min_y, point.y);
      *max_y = max(*max_y, point.y);
    }
  }
  fclose(points_file);
}

void DelaunayMesh::PlotPoints() const {
  double min_x;
  double max_x;
  double min_y;
  double max_y;
  SavePointsToFile("data_points", &min_x, &max_x, &min_y, &max_y);
  // Defines the plotting commands.
  FILE* plot_file = fopen("graph.conf", "w");
  fprintf(plot_file, "set xtic auto\nset ytic auto\n");
  fprintf(plot_file, "set title \"Delaunay\" \n set xlabel \"X\"\n");
  fprintf(plot_file, "set ylabel \"Y\" \nset pointsize 1\n");

  double xs = (max_x - min_x) * 0.2, ys = (max_y - min_y) * 0.2;

  fprintf(plot_file, "set xr [%lf : %lf]\n", min_x - xs, max_x + xs);
  fprintf(plot_file, "set yr [%lf : %lf]\n", min_y - ys, max_y + ys);
  fprintf(plot_file, "plot \"data_points\" using 1:2 with points pt 2"
                     "title \"Points\" ");
  fclose(plot_file);
  system("killall gnuplot");
  system("gnuplot -persist graph.conf");
}

bool DelaunayMesh::IsValidTriangle(const DelaunayTriangle& triangle) const {
  assert(triangle.triangle_child_ids.empty());
  bool contains_inf_points = false;
  for (const int point_id: triangle.point_ids) {
    contains_inf_points |= IsInfinitePoint(point_id);
  }
  if (!contains_inf_points) {
    for (int point_id = 0; point_id < _points.size(); ++point_id) {
      if (!IsInfinitePoint(point_id) &&
          IsInsideCircumcircle(_points[point_id],
                                GetPoint(triangle.point_ids[0]),
                                GetPoint(triangle.point_ids[1]),
                                GetPoint(triangle.point_ids[2]))) {
        return false;
      }
    }
  }
  return true;
}

bool DelaunayMesh::IsValidDelaunay() const {
  printf("Checking Delaunay Correctness ...\n");
  for (const DelaunayTriangle& triangle: _triangles) {
    // Only checks with triangle leaves.
    if (triangle.triangle_child_ids.empty()) {
      if(!IsValidTriangle(triangle)) {
        printf("\nCheck done: Delaunay is incorrect!!\n\n");
        return false;
      }
   }
  }
  printf("\nCheck done: Delaunay is correct!!\n\n");
  return true;
}


void DelaunayMesh::PlotTriangulation(bool wrong_circles,
                                     bool rand_circles,
                                     bool display_result) const {
  double min_x;
  double max_x;
  double min_y;
  double max_y;
  SavePointsToFile("data_points", &min_x, &max_x, &min_y, &max_y);

  FILE* plot_file = fopen("graph.conf", "w");

  fprintf(plot_file, "set size ratio -1\n"); // same units size
  fprintf(plot_file, "set xtic auto\nset ytic auto\n");
  fprintf(plot_file, "set title \"Delaunay\" \n set xlabel \"X\"\n");
  fprintf(plot_file, "set ylabel \"Y\" \nset pointsize 1\n");

  double xs = (max_x - min_x) * 0.1, ys = (max_y - min_y) * 0.1;

  fprintf(plot_file, "set xr [%lf : %lf]\n", min_x - xs, max_x + xs);
  fprintf(plot_file, "set yr [%lf : %lf]\n", min_y - ys, max_y + ys);
  fprintf(plot_file, "set style line 1 lc rgb \"#8D8D8D\"\n");
  fprintf(plot_file, "set style line 2 lc rgb \"#0000FF\"\n");
  fprintf(plot_file, "set style arrow 1 nohead ls 1\n");
  fprintf(plot_file, "set style arrow 2 nohead ls 2\n");

  // Contains the center of the circle and radius.
  FILE* circles_file = fopen("data_circles","w");
  // Contains the points that generates the circle.
  FILE* circle_points_file = fopen("data_circle_points","w");
  // Contains the data of invalid circles.
  FILE* invalid_circles_file = fopen("data_invalid_circles","w");
  // Contains the data of the points of invalid circles.
  FILE* invalid_circle_points_file = fopen("data_invalid_circle_points","w");
  // Checks if there exist invalid circles.
  bool invalid_circles = false;

  for (const DelaunayTriangle& triangle: _triangles) {
    // Plots only the leaves.
    if(triangle.triangle_child_ids.empty()) {
      // Plots connectivity among points.
      for (int i = 0; i < kTriangleSides; ++i) {
        int beg_point_id = triangle.point_ids[i];
        int end_point_id = triangle.point_ids[(i + 1) % 3];
        if (!IsInfinitePoint(beg_point_id) &&
            !IsInfinitePoint(end_point_id)) {
          const Point& beg_point = GetPoint(beg_point_id);
          const Point& end_point = GetPoint(end_point_id);
          if (beg_point == end_point) {
            fprintf(stderr, "%.6f %.6f\n%.6f %.6f\n\n", beg_point.x,
                    beg_point.y, end_point.x, end_point.y);
            fprintf(stderr, "%d %d\n", beg_point_id, end_point_id);
          }
          assert (beg_point != end_point);
          // Makes sure that the edge is printed only once.
          if (beg_point_id < end_point_id) {
            // Plot file is linear in the number of points.
            if (_constrains.count(minMaxPair(beg_point_id, end_point_id))) {
              fprintf(plot_file, "set arrow from %lf,%lf to %lf,%lf as 2\n",
                      beg_point.x, beg_point.y, end_point.x, end_point.y);
            } else {
              fprintf(plot_file, "set arrow from %lf,%lf to %lf,%lf as 1\n",
                  beg_point.x, beg_point.y, end_point.x, end_point.y);
            }
          }
        }
      }
      bool real_triangle = true;
      for (const int& point_id: triangle.point_ids) {
        real_triangle |= !IsInfinitePoint(point_id);
      }
      // Plots circumcircles.
      if (real_triangle) {
        // Plots 4% of the circumcircles.
        if(rand_circles && rand()%100 <= 4) {
          Point center;
          double radius;
          Circumcircle(GetPoint(triangle.point_ids[0]),
                       GetPoint(triangle.point_ids[1]),
                       GetPoint(triangle.point_ids[2]), &center, &radius);
          fprintf(circles_file, "%.4lf %.4lf %.4lf\n", center.x, center.y,
                  radius);

          for (const int& point_id: triangle.point_ids) {
            const Point& point = GetPoint(point_id);
            fprintf(circle_points_file, "%.4lf %.4lf\n", point.x, point.y);
          }
        }

        if(wrong_circles) {
          if (!IsValidTriangle(triangle)) {
            invalid_circles = true;
            Point center;
            double radius;
            Circumcircle(GetPoint(triangle.point_ids[0]),
                         GetPoint(triangle.point_ids[1]),
                         GetPoint(triangle.point_ids[2]), &center, &radius);
            fprintf(invalid_circles_file, "%.4lf %.4lf %.4lf\n", center.x,
                    center.y, radius);

            for (const int& point_id: triangle.point_ids) {
              const Point& point = GetPoint(point_id);
              fprintf(invalid_circle_points_file, "%.4lf %.4lf\n", point.x,
                      point.y);
            }
          }
        }
      }
    }
  }
  fclose(circles_file);
  fclose(circle_points_file);
  fclose(invalid_circles_file);
  fclose(invalid_circle_points_file);

  // Points plot command.
  fprintf(plot_file, "plot \"data_points\" using 1:2 with dots lc rgb \"black\""
         " title \"Points\"");

  if(invalid_circles) {
    // Invalid circles plot command.
    fprintf(plot_file, ", \\\n     \"data_invalid_circles\" with circles lw 1"
            " lc rgb \"red\", \\\n");

    // Invalid circle points plot command.
    fprintf(plot_file, "     \"data_invalid_circle_points\" with points pt 7"
            " lc rgb \"red\"");
  }

  if(rand_circles) {
    // Rand circles plot command.
    fprintf(plot_file, ", \\\n     \"data_circles\" with circles lw 1 lc rgb"
            " \"green\", \\\n");

    // Rand circle points plot command.
    fprintf(plot_file, "     \"data_circle_points\" with points pt 7 lc rgb"
            " \"green\"");
  }
  //fprintf(plot_file, "\n");
  fclose(plot_file);

  if (display_result) {
    system("killall gnuplot");
    system("gnuplot -persist graph.conf");
  }
}

void DelaunayMesh::GetSafeRegion(bool gnu_plot, int height, int width,
                                 Mat* safe_region) const {
  FILE* safe_circles_file;
  if (gnu_plot) {
    safe_circles_file = fopen("data_safe_circles", "w");
  }
  *safe_region = Mat::zeros(height, width, CV_8UC1);
  // Appends the safe circle information to the plot configuartion file.
  for (int triangle_id = 0; triangle_id < _triangles.size(); ++triangle_id) {
    const DelaunayTriangle& triangle = GetTriangleVal(triangle_id);
    // Plots only the leaves.
    if(triangle.triangle_child_ids.empty()) {
      // Plots connectivity among points.
      for (int side_id = 0; side_id < kTriangleSides; ++side_id) {
        int I = triangle.point_ids[side_id];
        int J = triangle.point_ids[(side_id + 1) % 3];
        int K = triangle.point_ids[(side_id + 2) % 3];
        if (!IsInfinitePoint(I) &&
            !IsInfinitePoint(J)) {
          const Point& beg_point = GetPoint(I);
          const Point& end_point = GetPoint(J);
          if (beg_point == end_point) {
            fprintf(stderr, "%.6f %.6f\n%.6f %.6f\n\n", beg_point.x,
                beg_point.y, end_point.x, end_point.y);
            fprintf(stderr, "%d %d\n", I, J);
          }
          assert (beg_point != end_point);
          // Makes sure that the edge is printed only once.
          if (I < J) {
            if (_constrains.count(minMaxPair(I, J))) {
              int neighbor_id = triangle.triangle_neighbor_ids[side_id];
              assert (GetTriangleVal(neighbor_id).triangle_child_ids.empty());
              int neighbor_side =  NeighborSide(neighbor_id, triangle_id);
              int L = GetTriangleVal(neighbor_id)
                      .point_ids[(neighbor_side + 2) %3];
              const Point pI = GetPoint(I);
              const Point pJ = GetPoint(J);
              const Point pK = GetPoint(K);
              const Point pL = GetPoint(L);
              Point center = (pI + pJ) / 2;
              double radius = center.Dist(pI);
              bool ok = true;
              if (!IsInfinitePoint(K)) {
                ok &= Cmp(center.Dist(pK), radius) >= 0;
              }
              if (ok && !IsInfinitePoint(L)) {
                ok &= Cmp(center.Dist(pL), radius) >= 0;
              }
              // Circle with diameter on the edge.
              if (ok) {
                if (gnu_plot) {
                  fprintf(safe_circles_file, "%.4lf %.4lf %.4lf\n", center.x,
                          center.y, radius);
                }
                Pixel center_pixel = center.ToPixel(height);
                // Draws the circle on safe_region.
                cv::circle(*safe_region,
                    cv::Point(center.x, height - 1 - center.y),
                    radius+2, cv::Scalar(kSafeRegion), -1);
                continue;
              }
              // Circles passing through circumcircles of adjacent neighbors.
              assert(!(IsInfinitePoint(K) && IsInfinitePoint(L)));
              radius = 1./0.;
              if (!IsInfinitePoint(K)) {
                Circumcircle(pI, pJ, pK, &center, &radius);
              }
              if (!IsInfinitePoint(L)) {
                Point qcenter;
                double qradius;
                Circumcircle(pI, pJ, pL, &qcenter, &qradius);
                if (Cmp(qradius, radius) < 0) {
                  radius = qradius;
                  center = qcenter;
                }
              }
              if (gnu_plot) {
                fprintf(safe_circles_file, "%.4lf %.4lf %.4lf\n", center.x,
                        center.y, radius);
              }
              Pixel center_pixel = center.ToPixel(height);
              // Draws the circle on safe_region.
              cv::circle(*safe_region,
                         cv::Point(center.x, height - 1 - center.y),
                         radius+2, cv::Scalar(kSafeRegion), -1);
            }
          }
        }
      }
    }
  }

  if (gnu_plot) {
    fclose(safe_circles_file);
    // Generates the plot configuration file for the triangle mesh.
    PlotTriangulation(false, false, false);

    // Adding the plot safe circles command.
    FILE* plot_file = fopen("graph.conf", "a");
    fprintf(plot_file, ", \\\n     \"data_safe_circles\" with circles lw 1"
        " lc rgb \"red\" fill solid noborder");
    fclose(plot_file);

    system("killall gnuplot");
    system("gnuplot -persist graph.conf");
  }
}

int DelaunayMesh::TriangleTreeSize() const {
  return  _triangles.size();
}

int DelaunayMesh::PointSetSize() const {
  return _points.size();
}

double DelaunayMesh::AverageInsertionCost() const {
  return ((double)_total_query_operations)/(_points.size()-3);
}

const cv::Mat& DelaunayMesh::GetSafeRegion() const {
  return _safe_region;
}

}  // namespace mesh_generation
