#pragma once

#include <cmath>
#include <cstddef>
#include <vector>

#include "common_geometry/math.hpp"
#include "common_geometry/types.hpp"

// Internal post-processing pass for the "marching squares, corner
// sharpened" extraction method. Not part of the public API.
namespace ns_ls2p::detail {

// Intersection of the infinite lines through (a,b) and (c,d), or false if
// they're parallel (within `kParallelEpsilon`).
inline bool LineIntersection(const ns_cg::Vec2d& a, const ns_cg::Vec2d& b,
                              const ns_cg::Vec2d& c, const ns_cg::Vec2d& d,
                              ns_cg::Vec2d& out) {
  const ns_cg::Vec2d dir1 = b - a;
  const ns_cg::Vec2d dir2 = d - c;
  const double denom = dir1.x() * dir2.y() - dir1.y() * dir2.x();
  constexpr double kParallelEpsilon = 1e-12;
  if (std::abs(denom) < kParallelEpsilon) return false;

  const double t1 = ((c.x() - b.x()) * dir2.y() - (c.y() - b.y()) * dir2.x()) / denom;
  out = b + t1 * dir1;
  return true;
}

// Post-processes a closed vertex loop: whenever two consecutive vertices
// B,C both have an interior angle within `angle_tolerance_degrees` of
// `target_angle_degrees` (135 degrees is marching squares' characteristic
// diagonal-chamfer angle for a 90-degree corner), replaces the pair with
// the single point where their neighboring edges (A->B and C->D extended)
// intersect -- recovering the sharp corner the chamfer cut off.
//
// This is a heuristic, not a general algorithm: it only fires on vertex
// pairs that look like a chamfer by this angle test, so it leaves loops
// with no such pairs (e.g. a smooth curve's vertices, which don't tend to
// land near a fixed angle) unchanged.
inline std::vector<ns_cg::Vec2d> SharpenCorners(
    const std::vector<ns_cg::Vec2d>& loop, double target_angle_degrees = 135.0,
    double angle_tolerance_degrees = 5.0) {
  const std::size_t n = loop.size();
  if (n < 4) return loop;

  std::vector<ns_cg::Vec2d> result;
  result.reserve(n);
  std::vector<bool> consumed(n, false);

  for (std::size_t i = 0; i < n; ++i) {
    if (consumed[i]) continue;
    const std::size_t b_idx = i;
    const std::size_t c_idx = (i + 1) % n;
    if (consumed[c_idx]) {
      result.push_back(loop[b_idx]);
      continue;
    }

    const ns_cg::Vec2d& a = loop[(b_idx + n - 1) % n];
    const ns_cg::Vec2d& b = loop[b_idx];
    const ns_cg::Vec2d& c = loop[c_idx];
    const ns_cg::Vec2d& d = loop[(c_idx + 1) % n];

    const double angle_b = ns_cg::InteriorAngleDegrees(a, b, c);
    const double angle_c = ns_cg::InteriorAngleDegrees(b, c, d);
    const bool looks_like_chamfer =
        std::abs(angle_b - target_angle_degrees) <= angle_tolerance_degrees &&
        std::abs(angle_c - target_angle_degrees) <= angle_tolerance_degrees;

    ns_cg::Vec2d sharp_corner;
    if (looks_like_chamfer && LineIntersection(a, b, c, d, sharp_corner)) {
      result.push_back(sharp_corner);
      consumed[c_idx] = true;  // B and C both collapse into sharp_corner
    } else {
      result.push_back(b);
    }
  }
  return result;
}

}  // namespace ns_ls2p::detail
