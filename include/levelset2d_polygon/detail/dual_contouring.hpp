#pragma once

#include <cmath>
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

// Internal dual contouring implementation. Not part of the public API.
//
// Unlike marching squares (which places a vertex on whichever grid edge the
// zero crossing falls on, then connects two such edge-points with a
// straight line through the cell), dual contouring places one vertex
// *inside* the cell, positioned via each crossing's estimated surface
// normal -- this is what lets it recover sharp features like a right-angle
// corner instead of always cutting a diagonal chamfer across it.
//
// This reuses marching squares' case table (same 16 cases, same saddle
// disambiguation) to decide *which* pairs of edge crossings form each
// boundary segment, and its already-validated topology (detail::
// LinkIntoLoops/SignedArea from marching_squares.hpp) to link segments
// into loops -- only the *position* assigned to each segment differs.
namespace ns_ls2p::detail {

// Central-difference gradient estimate at grid node (i,j) (one-sided at
// the boundary). Points in the direction of increasing field value, i.e.
// outward, since inside is negative.
inline ns_cg::Vec2d EstimateGradient(const ns_cg::Grid2d<double>& field,
                                      std::size_t i, std::size_t j) {
  double gx;
  if (i == 0) {
    gx = (field.at(1, j) - field.at(0, j)) / field.dx();
  } else if (i + 1 == field.nx()) {
    gx = (field.at(i, j) - field.at(i - 1, j)) / field.dx();
  } else {
    gx = (field.at(i + 1, j) - field.at(i - 1, j)) / (2.0 * field.dx());
  }

  double gy;
  if (j == 0) {
    gy = (field.at(i, 1) - field.at(i, 0)) / field.dy();
  } else if (j + 1 == field.ny()) {
    gy = (field.at(i, j) - field.at(i, j - 1)) / field.dy();
  } else {
    gy = (field.at(i, j + 1) - field.at(i, j - 1)) / (2.0 * field.dy());
  }
  return ns_cg::Vec2d(gx, gy);
}

inline std::vector<ns_cg::Vec2d> EstimateGradients(
    const ns_cg::Grid2d<double>& field) {
  const std::size_t nx = field.nx();
  const std::size_t ny = field.ny();
  std::vector<ns_cg::Vec2d> grads(nx * ny);
  for (std::size_t j = 0; j < ny; ++j)
    for (std::size_t i = 0; i < nx; ++i)
      grads[j * nx + i] = EstimateGradient(field, i, j);
  return grads;
}

// Interpolates the (unit) surface normal at the zero crossing along
// segment [a,b] with values va,vb and gradients ga,gb, using the same
// interpolation parameter as the crossing's position.
inline ns_cg::Vec2d InterpolateNormal(double va, const ns_cg::Vec2d& ga,
                                       double vb, const ns_cg::Vec2d& gb) {
  const double t = va / (va - vb);
  const ns_cg::Vec2d n = ga + t * (gb - ga);
  const double len = n.norm();
  if (len < 1e-12) return ns_cg::Vec2d(0.0, 0.0);  // degenerate; SolveQef falls back
  return n / len;
}

// Solves for the point x minimizing sum_k (n_k . (x - p_k))^2 for the two
// crossings (p0,n0), (p1,n1), via the 2x2 normal-equations solve. Falls
// back to the midpoint of p0,p1 if the system is near-singular (the two
// normals are close to parallel, i.e. the boundary is locally close to
// straight -- no sharp feature to resolve, and the QEF minimum is only
// weakly constrained along the tangent direction).
inline ns_cg::Vec2d SolveQef(const ns_cg::Vec2d& p0, const ns_cg::Vec2d& n0,
                             const ns_cg::Vec2d& p1, const ns_cg::Vec2d& n1) {
  const ns_cg::Mat2d A = n0 * n0.transpose() + n1 * n1.transpose();
  const ns_cg::Vec2d b = n0 * n0.dot(p0) + n1 * n1.dot(p1);

  const double det = A.determinant();
  if (std::abs(det) < 1e-9) return 0.5 * (p0 + p1);
  return A.inverse() * b;
}

// One boundary segment as both a marching-squares-style straight edge
// (reused purely for linking, via the already-validated topology) and the
// dual-contour vertex that should represent it in the final polygon.
struct DualSegment {
  ns_cg::Edge2d segment;
  ns_cg::Vec2d dual_vertex;
};

// Mirrors MarchCell's case table exactly, but returns dual-contour
// vertices (placed via SolveQef using each segment's 2 edge crossings and
// their interpolated normals) instead of straight-line segments.
inline std::vector<DualSegment> MarchCellDual(
    double v00, double v10, double v11, double v01, const ns_cg::Vec2d& p00,
    const ns_cg::Vec2d& p10, const ns_cg::Vec2d& p11, const ns_cg::Vec2d& p01,
    const ns_cg::Vec2d& g00, const ns_cg::Vec2d& g10, const ns_cg::Vec2d& g11,
    const ns_cg::Vec2d& g01) {
  const bool in00 = v00 < 0.0;
  const bool in10 = v10 < 0.0;
  const bool in11 = v11 < 0.0;
  const bool in01 = v01 < 0.0;
  const int case_index = (in00 ? 1 : 0) | (in10 ? 2 : 0) | (in11 ? 4 : 0) |
                          (in01 ? 8 : 0);
  if (case_index == 0 || case_index == 15) return {};

  struct Crossing {
    ns_cg::Vec2d point;
    ns_cg::Vec2d normal;
  };
  const auto B = [&] {
    return Crossing{Interpolate(p00, v00, p10, v10),
                     InterpolateNormal(v00, g00, v10, g10)};
  };
  const auto R = [&] {
    return Crossing{Interpolate(p10, v10, p11, v11),
                     InterpolateNormal(v10, g10, v11, g11)};
  };
  const auto T = [&] {
    return Crossing{Interpolate(p01, v01, p11, v11),
                     InterpolateNormal(v01, g01, v11, g11)};
  };
  const auto L = [&] {
    return Crossing{Interpolate(p00, v00, p01, v01),
                     InterpolateNormal(v00, g00, v01, g01)};
  };
  const auto Make = [](const Crossing& a, const Crossing& b) {
    return DualSegment{ns_cg::Edge2d(a.point, b.point),
                        SolveQef(a.point, a.normal, b.point, b.normal)};
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

inline std::vector<DualSegment> CollectSegmentsDual(
    const ns_cg::Grid2d<double>& field) {
  std::vector<DualSegment> segments;
  const std::size_t nx = field.nx();
  const std::size_t ny = field.ny();
  if (nx < 2 || ny < 2) return segments;

  const std::vector<ns_cg::Vec2d> grads = EstimateGradients(field);
  const auto grad = [&](std::size_t i, std::size_t j) { return grads[j * nx + i]; };

  for (std::size_t j = 0; j + 1 < ny; ++j) {
    for (std::size_t i = 0; i + 1 < nx; ++i) {
      const std::vector<DualSegment> cell_segments = MarchCellDual(
          field.at(i, j), field.at(i + 1, j), field.at(i + 1, j + 1),
          field.at(i, j + 1), field.position(i, j), field.position(i + 1, j),
          field.position(i + 1, j + 1), field.position(i, j + 1), grad(i, j),
          grad(i + 1, j), grad(i + 1, j + 1), grad(i, j + 1));
      segments.insert(segments.end(), cell_segments.begin(),
                       cell_segments.end());
    }
  }
  return segments;
}

inline std::vector<ns_cg::Polygon2d> ExtractPolygonsDualContouring(
    const ns_cg::Grid2d<double>& field) {
  const std::vector<DualSegment> dual_segments = CollectSegmentsDual(field);

  std::vector<ns_cg::Edge2d> segments;
  segments.reserve(dual_segments.size());
  std::map<std::pair<VertexKey, VertexKey>, ns_cg::Vec2d> dual_lookup;
  for (const auto& ds : dual_segments) {
    segments.push_back(ds.segment);
    dual_lookup[{KeyOf(ds.segment.GetStart()), KeyOf(ds.segment.GetEnd())}] =
        ds.dual_vertex;
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
      verts.push_back(dual_lookup.at(key));
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
