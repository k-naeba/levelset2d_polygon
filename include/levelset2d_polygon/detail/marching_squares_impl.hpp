#pragma once

#include <utility>
#include <vector>

#include "common_geometry/edge.hpp"
#include "common_geometry/grid.hpp"
#include "common_geometry/math.hpp"
#include "common_geometry/polygon.hpp"
#include "levelset2d_polygon/detail/marching_squares.hpp"

// Out-of-line definition of ns_ls2p::ExtractPolygons(), included from the
// bottom of marching_squares.hpp after the declaration is already visible.
namespace ns_ls2p {

inline std::vector<ns_cg::Polygon2d> ExtractPolygons(
    const ns_cg::Grid2d<double>& field) {
  const std::vector<ns_cg::Edge2d> segments = detail::CollectSegments(field);
  const std::vector<std::vector<ns_cg::Edge2d>> loops =
      detail::LinkIntoLoops(segments);

  std::vector<ns_cg::Polygon2d> result;
  std::vector<std::vector<ns_cg::Vec2d>> holes;

  for (const auto& loop : loops) {
    if (loop.empty()) continue;
    std::vector<ns_cg::Vec2d> verts;
    verts.reserve(loop.size());
    for (const auto& e : loop) verts.push_back(e.GetStart());

    if (detail::SignedArea(loop) > 0.0) {
      result.emplace_back(std::move(verts));
    } else {
      holes.push_back(std::move(verts));
    }
  }

  // Attach each hole to the outer polygon that contains it.
  for (auto& hole : holes) {
    if (hole.empty()) continue;
    for (auto& poly : result) {
      if (ns_cg::PointInPolygon(hole.front(), poly.GetOuter())) {
        poly.GetHoles().push_back(std::move(hole));
        break;
      }
    }
  }
  return result;
}

}  // namespace ns_ls2p
