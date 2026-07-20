#pragma once

#include <cstddef>
#include <map>
#include <tuple>
#include <vector>

#include "common_geometry/edge.hpp"
#include "common_geometry/grid.hpp"
#include "common_geometry/types.hpp"

// Internal marching-squares implementation used by ExtractPolygons(). Not
// part of the public API.
namespace ns_ls2p::detail {

// Linearly interpolates the zero-crossing point along segment [a,b], given
// scalar field values va, vb at those endpoints (must have opposite signs).
inline ns_cg::Vec2d Interpolate(const ns_cg::Vec2d& a, double va,
                                 const ns_cg::Vec2d& b, double vb) {
  const double t = va / (va - vb);
  return a + t * (b - a);
}

// Returns the 0-2 boundary segments for one grid cell, given corner values
// (negative = inside, matching common_geometry::SignedDistanceToPolygon's
// convention) and corner positions. Segments are oriented so inside is on
// the left of the direction of travel (same convention as
// rectilinear2d_boolean's boundary tracer): outer contours trace CCW,
// holes CW.
//
// The two "saddle" cases (opposite corners inside, adjacent corners
// outside, or vice versa) are ambiguous from the 4 corner samples alone;
// they're disambiguated using the average of the 4 corners as an
// approximation of the cell-center value.
inline std::vector<ns_cg::Edge2d> MarchCell(double v00, double v10, double v11,
                                             double v01, const ns_cg::Vec2d& p00,
                                             const ns_cg::Vec2d& p10,
                                             const ns_cg::Vec2d& p11,
                                             const ns_cg::Vec2d& p01) {
  const bool in00 = v00 < 0.0;
  const bool in10 = v10 < 0.0;
  const bool in11 = v11 < 0.0;
  const bool in01 = v01 < 0.0;
  const int case_index = (in00 ? 1 : 0) | (in10 ? 2 : 0) | (in11 ? 4 : 0) |
                          (in01 ? 8 : 0);
  if (case_index == 0 || case_index == 15) return {};

  const auto B = [&] { return Interpolate(p00, v00, p10, v10); };  // bottom
  const auto R = [&] { return Interpolate(p10, v10, p11, v11); };  // right
  const auto T = [&] { return Interpolate(p01, v01, p11, v11); };  // top
  const auto L = [&] { return Interpolate(p00, v00, p01, v01); };  // left

  switch (case_index) {
    case 1:
      return {ns_cg::Edge2d(B(), L())};
    case 2:
      return {ns_cg::Edge2d(R(), B())};
    case 3:
      return {ns_cg::Edge2d(R(), L())};
    case 4:
      return {ns_cg::Edge2d(T(), R())};
    case 6:
      return {ns_cg::Edge2d(T(), B())};
    case 7:
      return {ns_cg::Edge2d(T(), L())};
    case 8:
      return {ns_cg::Edge2d(L(), T())};
    case 9:
      return {ns_cg::Edge2d(B(), T())};
    case 11:
      return {ns_cg::Edge2d(R(), T())};
    case 12:
      return {ns_cg::Edge2d(L(), R())};
    case 13:
      return {ns_cg::Edge2d(B(), R())};
    case 14:
      return {ns_cg::Edge2d(L(), B())};
    case 5: {
      // Corners 00 and 11 inside; 10 and 01 outside.
      const double center = (v00 + v10 + v11 + v01) * 0.25;
      if (center < 0.0)
        return {ns_cg::Edge2d(B(), R()), ns_cg::Edge2d(T(), L())};
      return {ns_cg::Edge2d(B(), L()), ns_cg::Edge2d(T(), R())};
    }
    case 10: {
      // Corners 10 and 01 inside; 00 and 11 outside.
      const double center = (v00 + v10 + v11 + v01) * 0.25;
      if (center < 0.0)
        return {ns_cg::Edge2d(L(), B()), ns_cg::Edge2d(R(), T())};
      return {ns_cg::Edge2d(R(), B()), ns_cg::Edge2d(L(), T())};
    }
    default:
      return {};  // unreachable: case_index is a 4-bit value in [0,15]
  }
}

// Runs MarchCell over every cell of `field` and collects all segments.
inline std::vector<ns_cg::Edge2d> CollectSegments(
    const ns_cg::Grid2d<double>& field) {
  std::vector<ns_cg::Edge2d> segments;
  const std::size_t nx = field.nx();
  const std::size_t ny = field.ny();
  if (nx < 2 || ny < 2) return segments;

  for (std::size_t j = 0; j + 1 < ny; ++j) {
    for (std::size_t i = 0; i + 1 < nx; ++i) {
      const std::vector<ns_cg::Edge2d> cell_segments = MarchCell(
          field.at(i, j), field.at(i + 1, j), field.at(i + 1, j + 1),
          field.at(i, j + 1), field.position(i, j), field.position(i + 1, j),
          field.position(i + 1, j + 1), field.position(i, j + 1));
      segments.insert(segments.end(), cell_segments.begin(),
                       cell_segments.end());
    }
  }
  return segments;
}

// A hashable/orderable key for a Vec2d vertex, used to link segments by
// exact endpoint match. Segments from two cells sharing a grid edge always
// interpolate that edge's crossing from the same two corner values in the
// same order, so shared endpoints are bit-identical -- no epsilon needed.
struct VertexKey {
  double x;
  double y;
  friend bool operator<(const VertexKey& a, const VertexKey& b) {
    return std::tie(a.x, a.y) < std::tie(b.x, b.y);
  }
};

inline VertexKey KeyOf(const ns_cg::Vec2d& v) { return {v.x(), v.y()}; }

// Links a soup of segments into closed loops by following, from each
// vertex, its unique outgoing segment. Assumes every vertex has exactly
// one outgoing segment (i.e. no degenerate junction where multiple contour
// branches meet at a single point).
inline std::vector<std::vector<ns_cg::Edge2d>> LinkIntoLoops(
    const std::vector<ns_cg::Edge2d>& segments) {
  std::map<VertexKey, ns_cg::Edge2d> outgoing;
  for (const auto& s : segments) outgoing[KeyOf(s.GetStart())] = s;

  std::vector<std::vector<ns_cg::Edge2d>> loops;
  std::map<VertexKey, bool> visited;
  for (const auto& s : segments) {
    const VertexKey start_key = KeyOf(s.GetStart());
    if (visited[start_key]) continue;
    std::vector<ns_cg::Edge2d> loop;
    VertexKey v = start_key;
    while (!visited[v]) {
      visited[v] = true;
      const ns_cg::Edge2d& next = outgoing.at(v);
      loop.push_back(next);
      v = KeyOf(next.GetEnd());
    }
    loops.push_back(std::move(loop));
  }
  return loops;
}

// Shoelace signed area of a closed loop given as an ordered edge chain.
// Positive = counterclockwise.
inline double SignedArea(const std::vector<ns_cg::Edge2d>& loop) {
  double area = 0.0;
  for (const auto& e : loop) {
    area +=
        e.GetStart().x() * e.GetEnd().y() - e.GetEnd().x() * e.GetStart().y();
  }
  return 0.5 * area;
}

}  // namespace ns_ls2p::detail
