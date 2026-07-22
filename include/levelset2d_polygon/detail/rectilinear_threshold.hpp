#pragma once

#include <cstddef>
#include <map>
#include <tuple>
#include <vector>

#include "common_geometry/edge.hpp"
#include "common_geometry/grid.hpp"
#include "common_geometry/math.hpp"
#include "common_geometry/polygon.hpp"
#include "common_geometry/types.hpp"

// Internal implementation of the "rectilinear threshold" extraction method:
// binarize the level set onto the grid's own cells (no interpolation) and
// trace axis-aligned boundaries, the same technique rectilinear2d_boolean
// uses for its own (unrelated, BBox2d-set-derived) occupancy grids. Not
// part of the public API.
//
// Reimplemented locally (not depended upon from rectilinear2d_boolean's
// detail:: namespace, which is that project's own private implementation)
// -- and considerably simpler here, since Grid2d is already a uniform grid
// with no coordinate-compression step needed.
namespace ns_ls2p::detail {

struct ThresholdGridVertex {
  std::size_t col;
  std::size_t row;
  friend bool operator<(const ThresholdGridVertex& a,
                         const ThresholdGridVertex& b) {
    return std::tie(a.col, a.row) < std::tie(b.col, b.row);
  }
};

struct ThresholdGridEdge {
  ThresholdGridVertex from;
  ThresholdGridVertex to;
  std::size_t component;
};

inline constexpr std::size_t kNoComponent = static_cast<std::size_t>(-1);

// Cell (col,row) [corners at grid nodes (col,row)..(col+1,row+1)] is
// "inside" if the average of its 4 corner samples is negative.
inline bool CellInside(const ns_cg::Grid2d<double>& field, std::size_t col,
                        std::size_t row) {
  const double avg = (field.at(col, row) + field.at(col + 1, row) +
                       field.at(col + 1, row + 1) + field.at(col, row + 1)) *
                      0.25;
  return avg < 0.0;
}

inline bool IsFilled(const ns_cg::Grid2d<double>& field, std::size_t num_cols,
                      std::size_t num_rows, std::size_t col, std::size_t row) {
  if (col >= num_cols || row >= num_rows) return false;
  return CellInside(field, col, row);
}

// Labels each filled cell with the id of its 4-connected component
// (non-filled cells get kNoComponent). O(num_cells) flood fill.
inline std::vector<std::size_t> LabelComponents(const ns_cg::Grid2d<double>& field,
                                                  std::size_t num_cols,
                                                  std::size_t num_rows,
                                                  std::size_t& num_components) {
  std::vector<std::size_t> label(num_cols * num_rows, kNoComponent);
  num_components = 0;
  std::vector<std::size_t> stack;

  for (std::size_t start = 0; start < label.size(); ++start) {
    const std::size_t start_row = start / num_cols;
    const std::size_t start_col = start % num_cols;
    if (!CellInside(field, start_col, start_row) ||
        label[start] != kNoComponent)
      continue;
    label[start] = num_components;
    stack.push_back(start);
    while (!stack.empty()) {
      const std::size_t idx = stack.back();
      stack.pop_back();
      const std::size_t row = idx / num_cols;
      const std::size_t col = idx % num_cols;
      auto visit = [&](std::size_t r, std::size_t c) {
        if (!IsFilled(field, num_cols, num_rows, c, r)) return;
        const std::size_t nidx = r * num_cols + c;
        if (label[nidx] == kNoComponent) {
          label[nidx] = num_components;
          stack.push_back(nidx);
        }
      };
      if (row > 0) visit(row - 1, col);
      visit(row + 1, col);
      if (col > 0) visit(row, col - 1);
      visit(row, col + 1);
    }
    ++num_components;
  }
  return label;
}

inline std::vector<ThresholdGridEdge> CollectBoundaryEdges(
    const ns_cg::Grid2d<double>& field, std::size_t num_cols,
    std::size_t num_rows, const std::vector<std::size_t>& labels) {
  std::vector<ThresholdGridEdge> edges;
  for (std::size_t row = 0; row < num_rows; ++row) {
    for (std::size_t col = 0; col < num_cols; ++col) {
      if (!IsFilled(field, num_cols, num_rows, col, row)) continue;
      const std::size_t component = labels[row * num_cols + col];

      if (row == 0 || !IsFilled(field, num_cols, num_rows, col, row - 1))
        edges.push_back(ThresholdGridEdge{{col, row}, {col + 1, row}, component});
      if (!IsFilled(field, num_cols, num_rows, col + 1, row))
        edges.push_back(
            ThresholdGridEdge{{col + 1, row}, {col + 1, row + 1}, component});
      if (!IsFilled(field, num_cols, num_rows, col, row + 1))
        edges.push_back(
            ThresholdGridEdge{{col + 1, row + 1}, {col, row + 1}, component});
      if (col == 0 || !IsFilled(field, num_cols, num_rows, col - 1, row))
        edges.push_back(ThresholdGridEdge{{col, row + 1}, {col, row}, component});
    }
  }
  return edges;
}

inline std::vector<std::vector<ThresholdGridEdge>> LinkIntoLoops(
    const std::vector<ThresholdGridEdge>& edges) {
  std::map<ThresholdGridVertex, ThresholdGridEdge> outgoing;
  for (const auto& e : edges) outgoing[e.from] = e;

  std::vector<std::vector<ThresholdGridEdge>> loops;
  std::map<ThresholdGridVertex, bool> visited;
  for (const auto& e : edges) {
    if (visited[e.from]) continue;
    std::vector<ThresholdGridEdge> loop;
    ThresholdGridVertex v = e.from;
    while (!visited[v]) {
      visited[v] = true;
      const ThresholdGridEdge& next = outgoing.at(v);
      loop.push_back(next);
      v = next.to;
    }
    loops.push_back(std::move(loop));
  }
  return loops;
}

// Converts a loop's grid-index vertices to world points, merging
// consecutive collinear unit edges, and returns the resulting minimal
// ordered edge chain.
inline std::vector<ns_cg::Edge2d> SimplifyLoop(
    const std::vector<ThresholdGridEdge>& loop, const ns_cg::Grid2d<double>& field) {
  auto to_point = [&](const ThresholdGridVertex& v) {
    return field.position(v.col, v.row);
  };

  std::vector<ns_cg::Vec2d> raw;
  raw.reserve(loop.size());
  for (const auto& e : loop) raw.push_back(to_point(e.from));

  std::vector<ns_cg::Vec2d> corners;
  corners.reserve(raw.size());
  const std::size_t n = raw.size();
  for (std::size_t i = 0; i < n; ++i) {
    const ns_cg::Vec2d& prev = raw[(i + n - 1) % n];
    const ns_cg::Vec2d& cur = raw[i];
    const ns_cg::Vec2d& next = raw[(i + 1) % n];
    const ns_cg::Vec2d d1 = cur - prev;
    const ns_cg::Vec2d d2 = next - cur;
    if (d1.x() * d2.y() - d1.y() * d2.x() != 0.0 ||
        d1.x() * d2.x() + d1.y() * d2.y() < 0.0) {
      corners.push_back(cur);
    }
  }

  std::vector<ns_cg::Edge2d> edges;
  edges.reserve(corners.size());
  const std::size_t m = corners.size();
  for (std::size_t i = 0; i < m; ++i)
    edges.push_back(ns_cg::Edge2d(corners[i], corners[(i + 1) % m]));
  return edges;
}

inline double ThresholdSignedArea(const std::vector<ns_cg::Edge2d>& edges) {
  double area = 0.0;
  for (const auto& e : edges)
    area +=
        e.GetStart().x() * e.GetEnd().y() - e.GetEnd().x() * e.GetStart().y();
  return 0.5 * area;
}

inline std::vector<ns_cg::Polygon2d> ExtractPolygonsRectilinearThreshold(
    const ns_cg::Grid2d<double>& field) {
  std::vector<ns_cg::Polygon2d> result;
  if (field.nx() < 2 || field.ny() < 2) return result;

  const std::size_t num_cols = field.nx() - 1;
  const std::size_t num_rows = field.ny() - 1;

  std::size_t num_components = 0;
  const std::vector<std::size_t> labels =
      LabelComponents(field, num_cols, num_rows, num_components);

  const std::vector<ThresholdGridEdge> edges =
      CollectBoundaryEdges(field, num_cols, num_rows, labels);
  const std::vector<std::vector<ThresholdGridEdge>> loops = LinkIntoLoops(edges);

  result.resize(num_components);
  for (const auto& loop : loops) {
    const std::size_t component = loop.front().component;
    const std::vector<ns_cg::Edge2d> chain = SimplifyLoop(loop, field);

    std::vector<ns_cg::Vec2d> polygon;
    polygon.reserve(chain.size());
    for (const auto& e : chain) polygon.push_back(e.GetStart());

    if (ThresholdSignedArea(chain) > 0.0) {
      result[component].GetOuter() = std::move(polygon);
    } else {
      result[component].GetHoles().push_back(std::move(polygon));
    }
  }
  return result;
}

}  // namespace ns_ls2p::detail
