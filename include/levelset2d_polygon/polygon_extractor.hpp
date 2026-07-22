#pragma once

#include <memory>
#include <vector>

#include "common_geometry/grid.hpp"
#include "common_geometry/polygon.hpp"

namespace ns_ls2p {

// Selects which algorithm CreatePolygonExtractor() builds. See each
// concrete extractor's doc comment (detail/polygon_extractor_impl.hpp) for
// how it trades off accuracy, sharp-feature fidelity, and simplicity.
enum class ExtractionMethod {
  kMarchingSquares,                // baseline: see marching_squares.hpp
  kMarchingSquaresCornerSharpened, // baseline + a corner-sharpening post-process
  kDualContouring,                 // per-segment vertex placement using estimated normals
  kSurfaceNets,                    // per-segment vertex placement at the crossing midpoint
  kRectilinearThreshold,           // binarize + axis-aligned boundary tracing, no interpolation
};

// Common interface for algorithms that extract polygons from a 2D scalar
// field's zero level set (negative = inside, matching
// ns_cg::SignedDistanceToPolygon's convention).
class PolygonExtractor {
public:
  virtual ~PolygonExtractor() = default;
  virtual std::vector<ns_cg::Polygon2d> Extract(
      const ns_cg::Grid2d<double>& field) const = 0;
};

// Factory method: creates a PolygonExtractor implementing `method`.
std::unique_ptr<PolygonExtractor> CreatePolygonExtractor(ExtractionMethod method);

}  // namespace ns_ls2p

#include "levelset2d_polygon/detail/polygon_extractor_impl.hpp"
