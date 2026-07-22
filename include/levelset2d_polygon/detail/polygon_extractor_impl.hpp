#pragma once

#include <memory>
#include <stdexcept>
#include <vector>

#include "common_geometry/grid.hpp"
#include "common_geometry/polygon.hpp"
#include "levelset2d_polygon/detail/corner_sharpening.hpp"
#include "levelset2d_polygon/detail/dual_contouring.hpp"
#include "levelset2d_polygon/detail/marching_squares.hpp"
#include "levelset2d_polygon/detail/rectilinear_threshold.hpp"
#include "levelset2d_polygon/detail/surface_nets.hpp"
#include "levelset2d_polygon/marching_squares.hpp"

// Concrete PolygonExtractor implementations and the CreatePolygonExtractor
// factory, included from the bottom of polygon_extractor.hpp after
// PolygonExtractor/ExtractionMethod are already visible. Not part of the
// public API surface beyond what polygon_extractor.hpp re-declares.
namespace ns_ls2p::detail {

// Baseline: see marching_squares.hpp. Interpolates the zero crossing
// linearly along each cell edge; simple and fast, but a right-angle
// corner is always cut by a single diagonal segment (see
// analysis/corner_chamfer_analysis.cpp).
class MarchingSquaresExtractor : public PolygonExtractor {
public:
  std::vector<ns_cg::Polygon2d> Extract(
      const ns_cg::Grid2d<double>& field) const override {
    return ns_ls2p::ExtractPolygons(field);
  }
};

// Baseline marching squares, then a post-process (detail::SharpenCorners)
// that detects vertex pairs shaped like a diagonal chamfer (both near 135
// degrees) and collapses them back to the sharp corner their neighboring
// edges imply. Cheap to add on top of an existing marching squares
// pipeline, but it's a heuristic: it only acts on angle patterns that
// match a chamfer, so it won't help (or hurt) shapes without one.
class MarchingSquaresCornerSharpenedExtractor : public PolygonExtractor {
public:
  std::vector<ns_cg::Polygon2d> Extract(
      const ns_cg::Grid2d<double>& field) const override {
    std::vector<ns_cg::Polygon2d> polygons = ns_ls2p::ExtractPolygons(field);
    for (auto& polygon : polygons) {
      polygon.GetOuter() = SharpenCorners(polygon.GetOuter());
      for (auto& hole : polygon.GetHoles()) hole = SharpenCorners(hole);
    }
    return polygons;
  }
};

// Places each boundary vertex inside its cell (not on a cell edge) using
// the surface normal estimated from the field's gradient, via a small
// least-squares fit (QEF). Recovers sharp corners much better than plain
// marching squares since the vertex isn't constrained to a grid edge; the
// tradeoff is needing gradient estimates (here, finite differences on the
// sampled field) and a somewhat more involved algorithm.
class DualContouringExtractor : public PolygonExtractor {
public:
  std::vector<ns_cg::Polygon2d> Extract(
      const ns_cg::Grid2d<double>& field) const override {
    return ExtractPolygonsDualContouring(field);
  }
};

// A simpler sibling of DualContouringExtractor: each boundary segment's
// vertex is the midpoint of its two edge crossings, rather than a
// QEF-fitted point from estimated normals. No gradient estimation needed,
// simpler to reason about, but doesn't pull the vertex toward a sharp
// corner as effectively as the QEF solve does.
class SurfaceNetsExtractor : public PolygonExtractor {
public:
  std::vector<ns_cg::Polygon2d> Extract(
      const ns_cg::Grid2d<double>& field) const override {
    return ExtractPolygonsSurfaceNets(field);
  }
};

// Binarizes the field per grid cell (no interpolation at all) and traces
// axis-aligned cell boundaries, the same technique rectilinear2d_boolean
// uses for its own occupancy grids. Right angles are reproduced exactly
// (every edge is axis-aligned), but diagonal or curved boundaries come out
// blocky/staircased rather than smooth.
class RectilinearThresholdExtractor : public PolygonExtractor {
public:
  std::vector<ns_cg::Polygon2d> Extract(
      const ns_cg::Grid2d<double>& field) const override {
    return ExtractPolygonsRectilinearThreshold(field);
  }
};

}  // namespace ns_ls2p::detail

namespace ns_ls2p {

inline std::unique_ptr<PolygonExtractor> CreatePolygonExtractor(
    ExtractionMethod method) {
  switch (method) {
    case ExtractionMethod::kMarchingSquares:
      return std::make_unique<detail::MarchingSquaresExtractor>();
    case ExtractionMethod::kMarchingSquaresCornerSharpened:
      return std::make_unique<detail::MarchingSquaresCornerSharpenedExtractor>();
    case ExtractionMethod::kDualContouring:
      return std::make_unique<detail::DualContouringExtractor>();
    case ExtractionMethod::kSurfaceNets:
      return std::make_unique<detail::SurfaceNetsExtractor>();
    case ExtractionMethod::kRectilinearThreshold:
      return std::make_unique<detail::RectilinearThresholdExtractor>();
  }
  throw std::invalid_argument("CreatePolygonExtractor: unknown ExtractionMethod");
}

}  // namespace ns_ls2p
