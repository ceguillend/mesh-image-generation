#include "curve.h"
#include <list>
#include <map>
#include <queue>
#include <stack>
#include <utility>
#include <algorithm>

using cv::Mat;
using cv::Vec3b;
using cv::Scalar;
using std::vector;
using std::list;
using std::map;
using std::queue;
using std::stack;
using std::pair;

namespace mesh_generation {

namespace {
const int kEdge = 255;
const int kRowDir[] = {-1, 1, 0, 0, 1, 1, -1, -1};
const int kColDir[] = {0, 0, -1, 1, -1, 1, -1, 1};
const int kNumNeighbors = 8;
// Size of a plotted point in pixels. 
const int kPointSize = 3;

class Dsu {
 public:
  Dsu(){}

  Dsu(int size): dsu_(size){
    for (int i = 0; i < size; ++i) {
      dsu_[i] = i;
    }
  }

  int Find(int i) {
    return dsu_[i] = dsu_[i] == i ? i : Find(dsu_[i]);
  }

  void Join(int i, int j) {
    dsu_[Find(i)] = Find(j);
  }

 private:
  vector<int> dsu_;
};
}  // namespace

void FindPixelCurves(const Mat& canny_edge, int min_component_size,
                      int curve_length_threshold,
                      list< vector<Point> >* point_curves,
                      Mat* curves_plot) {

  vector<Pixel> vertices;
  for (int i = 0; i < canny_edge.rows; ++i) {
    for (int j = 0; j < canny_edge.cols; ++j) {
      if (canny_edge.at<uchar>(i, j) == kEdge) {
        vertices.push_back(Pixel(i, j));
      }
    }
  }
  sort(vertices.begin(), vertices.end());
  // MST: ortogonal_edges = 0, diagonal_edges = 1 
  Dsu dsu(vertices.size()); 
  // Tree.
  vector< vector<int> > tree(vertices.size());
  for (int edge_type = 0; edge_type < 2; ++edge_type) {
    for (int vertex_id = 0; vertex_id < vertices.size(); ++vertex_id) {
      for (int dir = 4 * edge_type, c = 4; c; ++dir, --c) {
        Pixel next (vertices[vertex_id].row + kRowDir[dir],
                      vertices[vertex_id].col + kColDir[dir]);
        vector<Pixel>::iterator it = lower_bound(vertices.begin(),
                                                  vertices.end(), next);
        if (it != vertices.end() && *it == next) {
          int next_id = it - vertices.begin();
          if (dsu.Find(vertex_id) != dsu.Find(next_id)) {
            dsu.Join(vertex_id, next_id);
            tree[vertex_id].push_back(next_id);
            tree[next_id].push_back(vertex_id);
          }
        }
      }
    }
  }
  // Plots MST.
  // Removes small components.
  vector<Vec3b> pixel_color(vertices.size());
  vector<bool> visited(vertices.size(), false);
  vector< list< vector<int> >::iterator > curve_ref(vertices.size());
  vector<bool> valid(vertices.size(), true);
  for (int vertex = 0; vertex < tree.size(); ++vertex) {
    if (!visited[vertex] && tree[vertex].size() == 1) {
      stack<int> dfs;
      dfs.push(vertex);
      int component_size = 0;
      list< vector<int> > curves;
      curve_ref[vertex] = curves.insert(curves.end(), vector<int>());
      pixel_color[vertex] = rand_color();
      while (!dfs.empty()) {
        int current = dfs.top();
        dfs.pop();
        ++component_size;
        visited[current] = true;
        curve_ref[current]->push_back(current);
        for (int next: tree[current]) {
          if (!visited[next]) {
            if (tree[current].size() > 2) { // splits
              pixel_color[next] = rand_color();
              curve_ref[next] = curves.insert(curves.end(), vector<int>());
            } else {
              pixel_color[next] = pixel_color[current];
              curve_ref[next] = curve_ref[current]; 
            }
            dfs.push(next);
          }
        }
      }
      if (component_size < min_component_size) {
        for (const vector<int>& curve: curves) {
          for (int vertex: curve) {
            valid[vertex] = false;
          }
        }
      } else {
        for (const vector<int>& curve: curves) {
          if (curve.size() < curve_length_threshold) {
            for (int vertex: curve) {
              valid[vertex] = false;
            }
          } else {
            point_curves->push_back(vector<Point>());
            for (int vertex: curve) {
              point_curves->back().push_back(vertices[vertex]
                                                  .ToPoint(canny_edge.rows));
            }
          }
        }
      }
    }
  }

  // Plots pixel curves
	*curves_plot = Mat::zeros(canny_edge.rows, canny_edge.cols, CV_8UC3);
  for (int id = 0; id < vertices.size(); ++id) {
    if (valid[id]) {
      curves_plot->at<Vec3b>(vertices[id].row,
                                vertices[id].col) = pixel_color[id];
    }
  }
}


void SimplifyPointCurves(const std::list<std::vector<Point>>& point_curves,
                        int image_height, int image_width,
                        int max_point_distance, int max_segment_length,
                        std::list< std::vector<Point> >* simplified_curves,
                        cv::Mat* simplified_curves_plot) {
  simplified_curves->clear();
  for (const vector<Point>& curve: point_curves) {
    if (curve.size() < 2) {
      continue;
    }
    vector<Point> simplified_curve;
    stack< pair<int, int> > segments;
    segments.push(pair<int, int>(0, curve.size()-1));

    while (!segments.empty()) {
      const int left = segments.top().first;
      const int right = segments.top().second;
      segments.pop();
      bool split = right - left > 1;
      if (split) {
        int split_pos;
        // find equation
        const double a = curve[right].y - curve[left].y;
        const double b = -(curve[right].x - curve[left].x);
        const double c = -(a * curve[left].x + b * curve[left].y);
        const double norm = sqrt(Sqr(a) + Sqr(b));
        double far_dist = -1;
        //find farthest point
        for (int pos = left + 1; pos < right; ++pos) {
          const double dist = fabs(a * curve[pos].x + b * curve[pos].y + c)
                                  / norm;	
          if (Cmp(far_dist, dist) < 0 ) {
            far_dist = dist;
            split_pos = pos;
          }
        }
        split = false;
        if (Cmp(max_point_distance, far_dist) < 0) {
          split = true;
        } else if (Cmp(max_segment_length,
                     curve[left].Dist(curve[right])) < 0) {
          split = true;
          split_pos = (left + right) / 2;
        }
        if (split) {
          segments.push(pair<int, int> (split_pos, right));
          segments.push(pair<int, int> (left, split_pos));
        }
      } 
      if (!split) {
        // Final segment.
        simplified_curve.push_back(curve[left]);
      }
    }
    simplified_curve.push_back(curve.back());
    simplified_curves->push_back(simplified_curve);
  }

  // Draws the simplified curves.
  *simplified_curves_plot = Mat::zeros(image_height , image_width, CV_8UC3);
  for (const vector<Point>& curve: *simplified_curves) {
    assert(curve.size() > 1);
    Vec3b rgb_color = rand_color();
		Scalar curve_color(rgb_color[0], rgb_color[1], rgb_color[2]);
    bool first = true;
    Pixel previous_pixel;
    for (const Point& point: curve) {
      Pixel pixel = point.ToPixel(image_height);
      circle(*simplified_curves_plot, cvPoint(pixel.col, pixel.row),
                kPointSize, curve_color);
      if (!first) {
        line(*simplified_curves_plot,
            cvPoint(previous_pixel.col, previous_pixel.row),
            cvPoint(pixel.col, pixel.row), curve_color);
      }
      previous_pixel = pixel;
      first = false;
    }
  }
}

}  // namespace mesh_generation
