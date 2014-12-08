#include "delaunay.h"
#include <assert.h>
#include <stdio.h>
#include <algorithm>
#include <stack>

using cv::Mat;
using cv::Vec3b;
using cv::Vec3d;
using std::list;
using std::max;
using std::make_pair;
using std::min;
using std::pair;
using std::stack;
using std::string;
using std::swap;
using std::unordered_map;
using std::vector;

const int kTriangleSides = 3;
const double kInfinite = 1./0.;
const int kP2Id = 0;
const int kP1Id = 1;
const int kNullId = -66666666;
const int kTriangleRoot = 0;
const bool kPlotSafeRegion = false;
const int kSafeRegion = 255;
/**
 * 4-pixel neighborhood.
 */
const vector<int> row_4nei =  {1,-1,0,0};
const vector<int> col_4nei =  {0,0,1,-1};
/**
 * Color neighborhood size.
 */
const int kColorNeighborhood = 3;

/**
 * Defines a Vec3b hashing.
 */
namespace std {
  template<>
    struct hash<Vec3b> {
      size_t operator ()(const Vec3b& color) const {
        int value = 0;
        for (int i = 0; i < 3; ++i) {
          value <<= 8;
          value += color.val[i];
        }
        return std::hash<int>()(value);
      }
    };
}  // namespace std

namespace {
/**
 * Canonical pair representation.
 */
pair<int, int> MinMaxPair(int a, int b) {
  return make_pair(min(a, b), max(a,b));
}

/**
 * Given an image and a set of pixels, finds the most frequent color. If
 * the set of pixels is empty returns false
 */
bool GetDominantColor(const Mat& source_image,
                      const vector<mesh_generation::Pixel>& pixel_group,
                      Vec3b* color) {
  if (pixel_group.empty()) {
    return false;
  }
  unordered_map<Vec3b, int> color_occurrences;
  for (const auto& pixel: pixel_group) {
    Vec3b pixel_color = source_image.at<Vec3b>(pixel.row, pixel.col);
    int var[3];
    for (var[0] = -kColorNeighborhood;
        var[0] <= kColorNeighborhood; ++var[0]) {
      for (var[1] = -kColorNeighborhood;
          var[1] <= kColorNeighborhood; ++var[1]) {
        for (var[2] = -kColorNeighborhood;
            var[2] <= kColorNeighborhood; ++var[2]) {
          int ncolor[3];
          bool valid = true;
          for (int i = 0; i < 3 && valid; ++i) {
            ncolor[i] = pixel_color.val[i] + var[i];
            if (ncolor[i] < 0 || ncolor[i] >= 256) {
              valid = false;
            }
          }
          if (valid) {
            ++color_occurrences[Vec3b(ncolor[0], ncolor[1], ncolor[2])];
          }
        }
      }
    }
  }
  int max_occurrence = 0;
  for (const auto& pixel: pixel_group) {
    Vec3b pixel_color = source_image.at<Vec3b>(pixel.row, pixel.col);
    const int& color_occurrence = color_occurrences[pixel_color];
    if (max_occurrence < color_occurrence) {
      max_occurrence = color_occurrence;
      *color = pixel_color;
    }
  }
  return true;
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

DelaunayMesh::DelaunayMesh() {}

DelaunayMesh::DelaunayMesh(int height, int width,
                           const list< vector<Point> >& simplified_curves,
                           int num_split_operations) {
  assert(num_split_operations >= 0);
  _total_query_operations = 0;
  Point lexicographically_bigger(width, height);
  // Creates the "outside-triangle" which is big enough to contain any triangle
  // inserted, the coordinates of the imaginary points are not actually used.
  int p2_id = PointId(Point(-kInfinite, kInfinite));
  int p1_id = PointId(Point(kInfinite, -kInfinite));
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
        _constrains[MinMaxPair(point1, point2)] = num_split_operations;
        _unsatisfied_constrains.insert(MinMaxPair(point1, point2));
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
    _constrains[MinMaxPair(new_point, it->first)] = _constrains[*it] - 1;
    _constrains[MinMaxPair(new_point, it->second)] = _constrains[*it] - 1;
    _constrains[*it] = 0;
    _unsatisfied_constrains.erase(it);
    //  Inserts the point to the delaunay triangulation.
    InsertPoint(new_point);
    ++num_inserted_points;
  }
  printf("Number of additional points inserted: %d\n", num_inserted_points);

  PairSet satisfied_constrains;
  GetSatisfiedConstrains(&satisfied_constrains);

  vector< pair<int, int> > unused_constrains;

  // Clears the unsatisfied constrains.
  for (const auto& constrain: _constrains) {
    if (satisfied_constrains.count(constrain.first) == 0) {
      unused_constrains.push_back(constrain.first);
    }
  }
  for (const auto& constrain: unused_constrains) {
    _constrains.erase(constrain);
  }
  GetSafeRegion(kPlotSafeRegion, height, width, &_safe_region);
}

DelaunayMesh::DelaunayMesh(int height, int width,
                           const vector< pair<Point,
                                              Point> >& segment_constrains) {
  _total_query_operations = 0;
  Point lexicographically_bigger(width, height);
  // Creates the "outside-triangle" which is big enough to contain any triangle
  // inserted, the coordinates of the imaginary points are not actually used.
  int p2_id = PointId(Point(-kInfinite, kInfinite));
  int p1_id = PointId(Point(kInfinite, -kInfinite));
  int p0_id = PointId(lexicographically_bigger);

  int outer_triangle_id = AllocateNewTriangle();
  GetTriangleRef(outer_triangle_id).point_ids = {p2_id, p0_id, p1_id};
  for (const pair<Point, Point>& constrain: segment_constrains) {
    int point1_id = PointId(constrain.first);
    int point2_id = PointId(constrain.second);
    // Does not allow splitting procedures.
    _constrains[MinMaxPair(point1_id, point2_id)] = 0;
  }
  vector<int> point_ids;
  for (int point_id = 3; point_id < _points.size(); ++point_id) {
    point_ids.push_back(point_id);
  }
  for (int i = 0; i < point_ids.size(); ++i) {
    swap(point_ids[i], point_ids[rand() % (i + 1)]);
  }
  for (int i = 0; i < point_ids.size(); ++i) {
    InsertPoint(point_ids[i]);
  }
  // Check constrains.
  PairSet satisfied_constrains;
  GetSatisfiedConstrains(&satisfied_constrains);
  int constrains_preserved = 0;
  for (const pair<Point, Point>& constrain: segment_constrains) {
    int point1_id = PointId(constrain.first);
    int point2_id = PointId(constrain.second);
    bool valid = satisfied_constrains.count(MinMaxPair(point1_id, point2_id));
    constrains_preserved += valid;
    /*
    assert(valid || !(std::cerr << width << " " << height << std::endl
                                << constrain.first.x << " "
                                << constrain.first.y << std::endl
                                << constrain.second.x << " "
                                << constrain.second.y << std::endl));
    */
  }
  // The assertion is not made to be constrains_preserved == segment_constrains.
  // since there can be double approximation errors (like cocircular points),
  // and it is possible that the point (width, height) might have not been
  // included in the set of original constrains thus breaking some other
  // constrains. The difference value is just an empirical estimation.
  assert(constrains_preserved <= segment_constrains.size());
  assert(segment_constrains.size() - constrains_preserved <= 10 ||
         !(std::cerr << constrains_preserved << " out of "
                     << segment_constrains.size() << std::endl));
  GetSafeRegion(kPlotSafeRegion, height, width, &_safe_region);
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

int DelaunayMesh::FindEnclosingTriangle(const Point& point) const {
  // Starts at the triangle root.
  int triangle_id = kTriangleRoot;

  while (!GetTriangleVal(triangle_id).triangle_child_ids.empty()) {
    int next_triangle_id = kNullId;
    for (int child_triangle_id :
              GetTriangleVal(triangle_id).triangle_child_ids) {
      //++_total_query_operations;
      if (IsInsideTriangle(child_triangle_id, point)) {
        next_triangle_id = child_triangle_id;
        break;
      }
    }
    assert(next_triangle_id != kNullId ||
           !(std::cerr << point.x << " " << point.y << std::endl));
    triangle_id = next_triangle_id;
  }
  return triangle_id;
}


void DelaunayMesh::UpdateNeighbor(int triangle_id, int old_neighbor,
                                  int new_neighbor) {
  if (triangle_id == kNullId) {
    return;
  }
  assert(GetTriangleVal(triangle_id).triangle_child_ids.empty());
  assert(GetTriangleVal(new_neighbor).triangle_child_ids.empty());
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
  assert(GetTriangleVal(triangle_id).triangle_child_ids.empty());
  int neighbor_id = GetTriangleVal(triangle_id)
    .triangle_neighbor_ids[side_index];
  if (neighbor_id == kNullId) {
    // It must be the outer triangle, does not need legalizations.
    bool outer = false;
    for (int point_id: GetTriangleVal(triangle_id).point_ids) {
      outer |= IsInfinitePoint(point_id);
    }
    assert(outer);
    return ;
  }
  assert(GetTriangleVal(neighbor_id).triangle_child_ids.empty());
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
  // K is inserted point, i.e. never infinite.
  assert(!IsInfinitePoint(K));
  // The case where both I and J are infinite can only happen in the,
  // outer triangle.
  assert(IsInfinitePoint(I) + IsInfinitePoint(J)  <= 1);

  bool ilegal_edge = false;
  if (!IsInfinitePoint(I) &&
      !IsInfinitePoint(J) &&
      !IsInfinitePoint(L)) {
    // All are normal points (counter-clockwise).
    ilegal_edge = IsInsideCircumcircle(GetPoint(L),
        GetPoint(I), GetPoint(J),GetPoint(K));
  } else {
    if(!IsInfinitePoint(L)) {
      // One of the elements of the edge is an infinite point.
      if(GetPoint(L) > GetPoint(K) ) {
        swap(L, K);
      }
      // GetPoint(L) < GetPoint(K), lexicografical orter.
      if (I > J) {
        // I: holds the infinite vertex.
        swap(I, J);
      }
      // Given that between K and L exists a segment IJ connecting to a
      // infinite point I that implies that J is between the lexicographycal
      // order of L and K,
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
      // L is a infinite point, it is always a legal edge.
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
    _unsatisfied_constrains.erase(MinMaxPair(K, L));
    const auto it = _constrains.find(MinMaxPair(I, J));
    if (it != _constrains.end() && it->second > 0) {
      if (!IsInfinitePoint(I) && !IsInfinitePoint(J)) {
        _unsatisfied_constrains.insert(MinMaxPair(I, J));
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
    const vector<int>& current_children = triangle_children[i];
    const vector<int>& oposite_children = triangle_children[i ^ 1];

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
      .erase(MinMaxPair(point_id,
            current_triangle.point_ids[current_side]));
    _unsatisfied_constrains
      .erase(MinMaxPair(point_id,
            current_triangle.point_ids[(current_side + 2) % 3]));
  }
  vector<int> edge_points = {
    GetTriangleVal(triangle_id).point_ids[side_index],
    GetTriangleVal(triangle_id).point_ids[(side_index + 1) % 3]
  };

  // Check removed edge, for the constrains processing.
  const pair<int, int> rm_edge = MinMaxPair(edge_points[0], edge_points[1]);
  const auto it = _constrains.find(rm_edge);
  if (it != _constrains.end()) {
    if (it->second > 0) {
      // replaces the constrain.
      for (int i = 0; i < edge_points.size(); ++i) {
        const pair<int, int> new_edge = MinMaxPair(point_id, edge_points[i]);
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
    // Update the neighbor's information and checks that the neighbor is a leaf.
    UpdateNeighbor(GetTriangleVal(triangle_id).triangle_neighbor_ids[i],
        triangle_id, triangle_child);
  }
  // Sets the triangle's children
  GetTriangleRef(triangle_id).triangle_child_ids = triangle_children;

  // Updates constrains.
  for (int i = 0; i < kTriangleSides; ++i) {
    _unsatisfied_constrains
      .erase(MinMaxPair(point_id, GetTriangleVal(triangle_id).point_ids[i]));
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
          const auto segment = MinMaxPair(beg_point_id, end_point_id);
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
                                     bool display_result,
                                     bool kill_other_plots) const {
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
            if (_constrains.count(MinMaxPair(beg_point_id, end_point_id))) {
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
    if (kill_other_plots) {
      system("killall gnuplot");
    }
    system("gnuplot -persist graph.conf");
  }
}

void DelaunayMesh::GetSafeCircle(int triangle_id, int side_id,
                                 Point* center, double* radius) const {
  assert (GetTriangleVal(triangle_id).triangle_child_ids.empty());
  const DelaunayTriangle& triangle = GetTriangleVal(triangle_id);
  int I = triangle.point_ids[side_id];
  int J = triangle.point_ids[(side_id + 1) % 3];
  int K = triangle.point_ids[(side_id + 2) % 3];
  assert(!IsInfinitePoint(I));
  assert(!IsInfinitePoint(J));
  int neighbor_id = triangle.triangle_neighbor_ids[side_id];
  assert (GetTriangleVal(neighbor_id).triangle_child_ids.empty() ||
      !(std::cerr << "Tr_id: " << neighbor_id << " / " << _triangles.size()
                  << " " << GetTriangleVal(neighbor_id).triangle_child_ids.size()
                  << " " << GetTriangleVal(neighbor_id).triangle_child_ids[0]
                  << " " << GetTriangleVal(neighbor_id).triangle_child_ids[1]
                  << " " << GetTriangleVal(neighbor_id).triangle_child_ids[2]
                  << std::endl));
  int neighbor_side =  NeighborSide(neighbor_id, triangle_id);
  int L = GetTriangleVal(neighbor_id).point_ids[(neighbor_side + 2) %3];
  const Point pI = GetPoint(I);
  const Point pJ = GetPoint(J);
  const Point pK = GetPoint(K);
  const Point pL = GetPoint(L);
  *center = (pI + pJ) / 2;
  *radius = center->Dist(pI);
  bool ok = true;
  if (!IsInfinitePoint(K)) {
    ok &= Cmp(center->Dist(pK), *radius) >= 0;
  }
  if (ok && !IsInfinitePoint(L)) {
    ok &= Cmp(center->Dist(pL), *radius) >= 0;
  }
  // Circle with diameter on the edge.
  if (ok) {
    return;
  }
  // Circles passing through circumcircles of adjacent neighbors.
  assert(!(IsInfinitePoint(K) && IsInfinitePoint(L)));
  *radius = 1./0.;
  if (!IsInfinitePoint(K)) {
    Circumcircle(pI, pJ, pK, center, radius);
  }
  if (!IsInfinitePoint(L)) {
    Point qcenter;
    double qradius;
    Circumcircle(pI, pJ, pL, &qcenter, &qradius);
    if (Cmp(qradius, *radius) < 0) {
      *radius = qradius;
      *center = qcenter;
    }
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
            if (_constrains.count(MinMaxPair(I, J))) {
              Point center;
              double radius;
              GetSafeCircle(triangle_id, side_id, &center, &radius);
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
    PlotTriangulation(false, false, false, false);

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

void DelaunayMesh::
GetConstrainedPoints(vector<Point>* constrained_points) const {
  vector<int> point_ids;
  for (const auto& constrain: _constrains) {
    point_ids.push_back(constrain.first.first);
    point_ids.push_back(constrain.first.second);
  }
  sort(point_ids.begin(), point_ids.end());
  point_ids.resize(unique(point_ids.begin(), point_ids.end())
      - point_ids.begin());
  constrained_points->clear();
  for (int point_id: point_ids) {
    constrained_points->push_back(GetPoint(point_id));
  }
}

void DelaunayMesh::GetRandSafePoints(double fraction,
    vector<Point>* safe_points) const {
  // Obtains the safe region
  int height = _safe_region.rows;
  int width = _safe_region.cols;
  vector< vector<bool> > visited(height, vector<bool> (width));

  for (int row = 0; row < height; ++row) {
    for (int col = 0; col < width; ++col) {
      visited[row][col] = _safe_region.at<uchar>(row, col) == kSafeRegion;
    }
  }

  safe_points->clear();
  for (int row = 0; row < height; ++row) {
    for (int col = 0; col < width; ++col) {
      if (!visited[row][col]) {
        vector<Point> inner_region_points;
        stack< pair<int, int> > dfs;
        dfs.push(make_pair(row, col));
        visited[row][col] = true;
        while (!dfs.empty()) {
          int qrow = dfs.top().first;
          int qcol = dfs.top().second;
          dfs.pop();
          inner_region_points.push_back(Pixel(qrow, qcol).ToPoint(height));
          for (int dir = 0; dir < row_4nei.size(); ++dir) {
            int nrow = qrow + row_4nei[dir];
            int ncol = qcol + col_4nei[dir];
            if (0 <= nrow && nrow < height &&
                0 <= ncol && ncol < width && !visited[nrow][ncol]) {
              visited[nrow][ncol] = true;
              dfs.push(make_pair(nrow, ncol));
            }
          }
        }
        int sample = round(fraction * inner_region_points.size());
        for (int i = 0; i < inner_region_points.size(); ++i) {
          int p = rand() % (i + 1);
          if (p < sample) {
            swap(inner_region_points[p], inner_region_points[i]);
          }
        }
        for (int i = 0; i < sample; ++i) {
          safe_points->push_back(inner_region_points[i]);
        }
      }
    }
  }
}


void DelaunayMesh::
     GetUnconstrainedAdjustedPoints(vector<Point>* adjusted_points) const {
  adjusted_points->clear();
  const int height = _safe_region.rows;
  const int width = _safe_region.cols;
  vector<bool> processed(_points.size(), false);
  for (int triangle_id = 0 ; triangle_id < _triangles.size(); ++triangle_id) {
    const DelaunayTriangle& triangle = GetTriangleVal(triangle_id);
    if(triangle.triangle_child_ids.empty()) {
      for (int side_id = 0; side_id < kTriangleSides; ++side_id) {
        const int point_id = triangle.point_ids[side_id];
        if ( movable_point_ids.count(point_id) == 1 && !processed[point_id]) {
          processed[point_id] = true;
          bool movable = true;
          // Generates the polygon enclosing the point in counter-clockwise
          // order.
          // Stores (triangle_id, side_id) of the current point.
          vector< pair<int, int> > neighbor_triangles;
          int qtriangle_id = triangle_id;
          int qside_id = side_id;
          do {
            neighbor_triangles.push_back(make_pair(qtriangle_id, qside_id));
            // Determines if the point is movable.
            int npoint_id = GetTriangleVal(qtriangle_id)
              .point_ids[(qside_id + 1) % 3];
            movable &= !IsInfinitePoint(npoint_id);
            movable &= _constrains.count(MinMaxPair(point_id, npoint_id)) == 0;
            // Move to the next triangle.
            qtriangle_id = GetTriangleVal(qtriangle_id)
              .triangle_neighbor_ids[qside_id];
            if (qtriangle_id == kNullId) {
              assert(!movable);
              break;
            }
            assert((qtriangle_id >= 0 && qtriangle_id < _triangles.size()) ||
                !(std::cerr  << "Tr_id: " << qtriangle_id << std::endl));
            assert(GetTriangleVal(qtriangle_id).triangle_child_ids.empty());
            qside_id = 0;
            while (GetTriangleVal(qtriangle_id).point_ids[qside_id]
                      != point_id) {
              ++qside_id;
              assert(qside_id < 3);
            }
          } while(movable && qtriangle_id != triangle_id);

          assert(movable);
          assert(neighbor_triangles.size() > 2);

          // Generates the actual points of the polygons.
          vector<Point> polygon;
          for (const pair<int, int>& neighbor_data: neighbor_triangles) {
            const vector<int>& point_ids = GetTriangleVal(neighbor_data.first)
                                           .point_ids;
            for (int i = 0; i < point_ids.size(); ++i) {
              assert(!IsInfinitePoint(point_ids[i]));
            }
            Point center;
            double radius;
            Circumcircle(GetPoint(point_ids[0]), GetPoint(point_ids[1]),
                GetPoint(point_ids[2]), &center, &radius);
            // The circumcircle can be outside the image region, but the polygon
            // will be clipped to fit inside.
            polygon.push_back(center);
          }
          // The set of line clips: (unit vector direction, line point), the
          // region of the polygon to the right of the vector direction will be
          // keep.
          vector < pair<Point, Point> > clips = {
            make_pair(Point(0, 1), Point(0, 0)),
            make_pair(Point(1, 0), Point(0, height-1)),
            make_pair(Point(0, -1), Point(width-1, height-1)),
            make_pair(Point(-1, 0), Point(width-1, 0))
          };
          // Clips the voronoi polygon by the constrained edges.
          for (const pair<int, int>& neighbor_data: neighbor_triangles) {
            const int ntriangle_id = neighbor_data.first;
            const int nside_id = neighbor_data.second;
            assert(GetTriangleVal(ntriangle_id).point_ids[nside_id]
                    == point_id);
            const int point2_id = GetTriangleVal(ntriangle_id)
                                  .point_ids[(nside_id + 1) % 3];
            const int point3_id = GetTriangleVal(ntriangle_id)
                                  .point_ids[(nside_id + 2) % 3];
            if (_constrains.count(MinMaxPair(point2_id, point3_id))) {
              // Constrains clips.
              Point safe_center;
              double safe_radius;
              GetSafeCircle(ntriangle_id, (nside_id + 1) % 3, &safe_center,
                            &safe_radius);
              const Point point2 = GetPoint(point2_id);
              const Point point3 = GetPoint(point3_id);
              const Point dir_vec = (point3 - point2).Unit();
              // Offsets double error precision.
              const Point clip_point = safe_center +
                                       dir_vec.Ort()*-(safe_radius + 0.01);
              clips.push_back(make_pair(dir_vec, clip_point));
            }
          }
          for (const auto& clip: clips) {
            const Point dir_vec = clip.first;
            const Point clip_point = clip.second;
            vector<Point> clipped_polygon;
            vector<bool> inside(polygon.size());
            for (int i = 0; i < polygon.size(); ++i) {
              inside[i] = Cmp(dir_vec ^ (polygon[i] - clip_point), 0) < 0;
            }
            for (int i = 0; i < polygon.size(); ++i) {
              if (inside[i]) {
                clipped_polygon.push_back(polygon[i]);
              }
              int j = (i + 1) % polygon.size();
              if (inside[i] ^ inside[j]) {
                Point v = polygon[j] - polygon[i];
                double t = ((polygon[i] - clip_point) * v.Ort()) /
                           (dir_vec * v.Ort());
                clipped_polygon.push_back(clip_point + dir_vec * t);
              }
            }

            assert(clipped_polygon.size() != 1  &&
                   clipped_polygon.size() != 2);

            polygon = clipped_polygon;
            if (polygon.empty()) {
              break;
            }
          }
          if (polygon.empty()) {
            adjusted_points->push_back(GetPoint(point_id));
          } else {
            // Checks coordinates values.
            for (Point point: polygon) {
              assert((Cmp(point.x, 0) >= 0 &&
                      Cmp(point.x, width-1) <= 0 &&
                      Cmp(point.y, 0) >= 0 &&
                      Cmp(point.y, height-1) <= 0) ||
                     !(std::cerr<< point.x << " " << point.y << std::endl));
            }
            adjusted_points->push_back(Centroid(polygon));
          }
        }
      }
    }
  }
}

void  DelaunayMesh::InsertSafePoints(const vector<Point>& safe_points) {
  int height = _safe_region.rows;
  int width = _safe_region.cols;
  // The constrains should be preserved after insertion of the safe points.
  PairSet initial_constrains;
  PairSet final_constrains;
  GetSatisfiedConstrains(&initial_constrains);
  for (const Point& point: safe_points) {
    assert((Cmp(point.x, 0) >= 0 &&
          Cmp(point.x, width-1) <= 0 &&
          Cmp(point.y, 0) >= 0 &&
          Cmp(point.y, height-1) <= 0) ||
        !(std::cerr<< point.x << " " << point.y << std::endl));
    int point_id = PointId(point);
    InsertPoint(point_id);
    movable_point_ids.insert(point_id);
  }
  // Doesn't count the precise number of constrains since it might be the case
  // that cocircular constrains inserted in different order might lead to
  // unexpected behaviors.
  GetSatisfiedConstrains(&final_constrains);
  int count = 0;
  for (const pair<int, int>& constrain: initial_constrains) {
    //assert(final_constrains.count(constrain) > 0);
    count += final_constrains.count(constrain) > 0;
  }
  // The value is defined experimentally.
  assert(initial_constrains.size() - count <= 10);
  printf("Number of inner points inserted: %d\n", (int)safe_points.size());
  printf("Constrains: %d / %d\n", count, (int)initial_constrains.size());
}

void DelaunayMesh::GetMinimalConstrainSet(
                   vector< pair<Point, Point> >* minimal_constrain_set) const {
  minimal_constrain_set->clear();
  for (const pair<pair<int, int>, int>& constrain: _constrains) {
    minimal_constrain_set->push_back(
                                make_pair(GetPoint(constrain.first.first),
                                          GetPoint(constrain.first.second)));
  }
}

void DelaunayMesh::GetMeshSolidInterpolation(const cv::Mat& source_image,
                                        cv::Mat* interpolation) const {
  printf("Generating the interpolation ...\n");
  int height = _safe_region.rows;
  int width = _safe_region.cols;
  assert(source_image.rows == height && source_image.cols == width);
  *interpolation = Mat::zeros(height, width, CV_8UC3);
  unordered_map<int, vector<Pixel> > pixel_partitions;

  // Groups the pixels in triangles.
  for (int row = 0; row < height; ++row) {
    for (int col = 0; col < width; ++col) {
      Pixel pixel = Pixel(row, col);
      int triangle_id = FindEnclosingTriangle(pixel.ToPoint(height));
      pixel_partitions[triangle_id].push_back(pixel);
    }
  }

  // For each triangle calculates the interpolation.
  for (const auto& triangle_partition: pixel_partitions) {
    const vector<Pixel>& pixels = triangle_partition.second;
    Vec3b color;
    assert(GetDominantColor(source_image, pixels, &color));
    for (const auto& pixel: pixels) {
      interpolation->at<Vec3b>(pixel.row, pixel.col) = color;
    }
  }
  printf("Interpolation completed!!\n\n");
}

void DelaunayMesh::GetMeshLinearInterpolation(const cv::Mat& source_image,
    cv::Mat* interpolation) const {
  printf("Generating the interpolation ...\n");
  int height = _safe_region.rows;
  int width = _safe_region.cols;
  assert(source_image.rows == height && source_image.cols == width);
  *interpolation = Mat::zeros(height, width, CV_8UC3);
  unordered_map<int, vector<Pixel> > pixel_partitions;

  // Groups the pixels in triangles.
  for (int row = 0; row < height; ++row) {
    for (int col = 0; col < width; ++col) {
      Pixel pixel = Pixel(row, col);
      int triangle_id = FindEnclosingTriangle(pixel.ToPoint(height));
      pixel_partitions[triangle_id].push_back(pixel);
    }
  }

  // For each triangle calculates the interpolation.
  for (const auto& triangle_partition: pixel_partitions) {
    const int triangle_id  = triangle_partition.first;
    const vector<Pixel>& pixels = triangle_partition.second;
    const vector<int>& point_ids = GetTriangleVal(triangle_id).point_ids;

    // Vertices properties.
    vector<Point> points;
    for (int i = 0; i < point_ids.size(); ++i) {
      if (!IsInfinitePoint(point_ids[i])) {
        points.push_back(GetPoint(point_ids[i]));
      }
    }
    vector< vector<Pixel> > pixel_groups(points.size(), vector<Pixel>());
    for (const Pixel& pixel: pixels) {
      int closest_point = -1;
      double distance = 1./0.;
      for (int i = 0; i < points.size(); ++i) {
        double qdistance = points[i].Dist(pixel.ToPoint(height));
        if (qdistance < distance) {
          distance = qdistance;
          closest_point = i;
        }
      }
      assert(closest_point != -1);
      pixel_groups[closest_point].push_back(pixel);
    }
    vector<Vec3b> colors;
    // Determines the color of each triangle vertex.
    for (int i = 0; i < points.size(); ++i ) {
      Vec3b color;
      if (GetDominantColor(source_image, pixel_groups[i], &color)) {
        colors.push_back(color);
      }
    }
    while (colors.size() != points.size()) {
      colors.push_back(colors.back());
    }
    // Assigns the colors to each pixel.
    for (const Pixel& pixel: pixels) {
      Point point = pixel.ToPoint(height);
      if (points.size() == 1) {
        interpolation->at<Vec3b>(pixel.row, pixel.col) = colors[0];
      } else if (points.size() == 2) {
        double d1 = point.Dist(points[0]);
        double d2 = point.Dist(points[1]);
        double f1 = d1 / (d1 + d2);
        double f2 = d2 / (d1 + d2);
        interpolation->at<Vec3b>(pixel.row, pixel.col) = Vec3d(colors[0]) * f1 +
                                                         Vec3d(colors[1]) * f2;
      } else if (points.size() == 3) {
        // Barycentric coordinates of the point.
        double bc[3];
        Barycentric(point, points[0], points[1], points[2], bc, bc + 1,
                    bc + 2);
        interpolation->at<Vec3b>(pixel.row, pixel.col) =
            Vec3d(colors[0]) * bc[0] +
            Vec3d(colors[1]) * bc[1] +
            Vec3d(colors[2]) * bc[2];
      } else {
        assert (false);
      }
    }
  }
  printf("Interpolation completed!!\n\n");
}
}  // namespace mesh_generation
