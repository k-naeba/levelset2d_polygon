#pragma once

#include <cstddef>
#include <map>
#include <utility>
#include <vector>

#include "common_geometry/edge.hpp"
#include "common_geometry/grid.hpp"
#include "common_geometry/math.hpp"
#include "common_geometry/polygon.hpp"
#include "common_geometry/types.hpp"
#include "levelset2d_polygon/detail/marching_squares.hpp"

// Internal surface nets implementation. Not part of the public API.
//
// A simpler sibling of dual contouring: rather than solving a QEF from
// estimated surface normals, each boundary segment's vertex is simply the
// midpoint of its two edge crossings -- no gradient estimation needed at
// all. Reuses marching squares' case table for segment topology (same 16
// cases, same saddle disambiguation) and its already-validated
// LinkIntoLoops/SignedArea for linking, exactly as dual contouring does;
// only the vertex-placement rule differs.
namespace ns_ls2p::detail {

// One boundary segment as both a marching-squares-style straight edge
// (reused purely for linking) and the surface-nets vertex that should
// represent it in the final polygon.
struct SurfaceNetsSegment {
  ns_cg::Edge2d segment;
  ns_cg::Vec2d vertex;
};

// Mirrors MarchCell's case table exactly, but returns the midpoint of each
// segment's 2 edge crossings instead of the straight-line segment itself.
inline std::vector<SurfaceNetsSegment> MarchCellSurfaceNets(
    double v00, double v10, double v11, double v01, const ns_cg::Vec2d& p00,
    const ns_cg::Vec2d& p10, const ns_cg::Vec2d& p11, const ns_cg::Vec2d& p01) {
  const bool in00 = v00 < 0.0;
  const bool in10 = v10 < 0.0;
  const bool in11 = v11 < 0.0;
  const bool in01 = v01 < 0.0;
  const int case_index = (in00 ? 1 : 0) | (in10 ? 2 : 0) | (in11 ? 4 : 0) |
                          (in01 ? 8 : 0);
  if (case_index == 0 || case_index == 15) return {};

  const auto B = [&] { return Interpolate(p00, v00, p10, v10); };
  const auto R = [&] { return Interpolate(p10, v10, p11, v11); };
  const auto T = [&] { return Interpolate(p01, v01, p11, v11); };
  const auto L = [&] { return Interpolate(p00, v00, p01, v01); };
  const auto Make = [](const ns_cg::Vec2d& a, const ns_cg::Vec2d& b) {
    return SurfaceNetsSegment{ns_cg::Edge2d(a, b), 0.5 * (a + b)};
  };

  switch (case_index) {
    case 1:
      return {Make(B(), L())};
    case 2:
      return {Make(R(), B())};
    case 3:
      return {Make(R(), L())};
    case 4:
      return {Make(T(), R())};
    case 6:
      return {Make(T(), B())};
    case 7:
      return {Make(T(), L())};
    case 8:
      return {Make(L(), T())};
    case 9:
      return {Make(B(), T())};
    case 11:
      return {Make(R(), T())};
    case 12:
      return {Make(L(), R())};
    case 13:
      return {Make(B(), R())};
    case 14:
      return {Make(L(), B())};
    case 5: {
      const double center = (v00 + v10 + v11 + v01) * 0.25;
      if (center < 0.0) return {Make(B(), R()), Make(T(), L())};
      return {Make(B(), L()), Make(T(), R())};
    }
    case 10: {
      const double center = (v00 + v10 + v11 + v01) * 0.25;
      if (center < 0.0) return {Make(L(), B()), Make(R(), T())};
      return {Make(R(), B()), Make(L(), T())};
    }
    default:
      return {};  // unreachable: case_index is a 4-bit value in [0,15]
  }
}

inline std::vector<SurfaceNetsSegment> CollectSegmentsSurfaceNets(
    const ns_cg::Grid2d<double>& field) {
  std::vector<SurfaceNetsSegment> segments;
  const std::size_t nx = field.nx();
  const std::size_t ny = field.ny();
  if (nx < 2 || ny < 2) return segments;

  for (std::size_t j = 0; j + 1 < ny; ++j) {
    for (std::size_t i = 0; i + 1 < nx; ++i) {
      const std::vector<SurfaceNetsSegment> cell_segments = MarchCellSurfaceNets(
          field.at(i, j), field.at(i + 1, j), field.at(i + 1, j + 1),
          field.at(i, j + 1), field.position(i, j), field.position(i + 1, j),
          field.position(i + 1, j + 1), field.position(i, j + 1));
      segments.insert(segments.end(), cell_segments.begin(),
                       cell_segments.end());
    }
  }
  return segments;
}

inline std::vector<ns_cg::Polygon2d> ExtractPolygonsSurfaceNets(
    const ns_cg::Grid2d<double>& field) {
  const std::vector<SurfaceNetsSegment> sn_segments = CollectSegmentsSurfaceNets(field);

  std::vector<ns_cg::Edge2d> segments;
  segments.reserve(sn_segments.size());
  std::map<std::pair<VertexKey, VertexKey>, ns_cg::Vec2d> vertex_lookup;
  for (const auto& s : sn_segments) {
    segments.push_back(s.segment);
    vertex_lookup[{KeyOf(s.segment.GetStart()), KeyOf(s.segment.GetEnd())}] =
        s.vertex;
  }

  // Reuse marching squares' already-validated linking/orientation logic --
  // only the final vertex positions differ.
  const std::vector<std::vector<ns_cg::Edge2d>> loops = LinkIntoLoops(segments);

  std::vector<ns_cg::Polygon2d> result;
  std::vector<std::vector<ns_cg::Vec2d>> holes;

  for (const auto& loop : loops) {
    if (loop.empty()) continue;
    std::vector<ns_cg::Vec2d> verts;
    verts.reserve(loop.size());
    for (const auto& e : loop) {
      const auto key = std::make_pair(KeyOf(e.GetStart()), KeyOf(e.GetEnd()));
      verts.push_back(vertex_lookup.at(key));
    }

    if (SignedArea(loop) > 0.0) {
      result.emplace_back(std::move(verts));
    } else {
      holes.push_back(std::move(verts));
    }
  }

  ns_cg::AttachHolesByContainment(result, std::move(holes));
  return result;
}

}  // namespace ns_ls2p::detail
