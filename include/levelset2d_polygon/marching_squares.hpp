#pragma once

#include <vector>

#include "common_geometry/grid.hpp"
#include "common_geometry/polygon.hpp"

namespace ns_ls2p {

// Extracts the zero level set of `field` (e.g. as produced by
// ns_cg::BuildLevelSet) via marching squares, as one ns_cg::Polygon2d per
// connected component: its outer contour, plus the contour of every hole
// it encloses. Uses the same sign convention as
// ns_cg::SignedDistanceToPolygon: negative = inside.
//
// Limitation: components that touch another component (or themselves) at a
// single point are not supported and may produce an incorrect trace (same
// caveat as rectilinear2d_boolean's ExtractContours()).
std::vector<ns_cg::Polygon2d> ExtractPolygons(const ns_cg::Grid2d<double>& field);

}  // namespace ns_ls2p

#include "levelset2d_polygon/detail/marching_squares_impl.hpp"
