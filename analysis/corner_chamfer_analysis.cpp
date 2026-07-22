// Analysis (not a usage example): compares all 4 polygon extraction
// methods (see polygon_extractor.hpp) on how well they reproduce a sharp
// right-angle corner. Baseline marching squares always replaces a 90
// degree corner with a diagonal cut -- two new vertices near the true
// corner, each at a 135-degree angle, regardless of grid resolution (only
// the cut's size shrinks). This program measures the same corner across
// the other 3 methods too, at a coarse and a fine resolution, and exports
// SVGs comparing them directly.

#include <algorithm>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "common_geometry/levelset.hpp"
#include "common_geometry/math.hpp"
#include "common_geometry/svg.hpp"
#include "levelset2d_polygon/levelset2d_polygon.hpp"

using ns_cg::BuildLevelSet;
using ns_cg::Grid2d;
using ns_cg::InteriorAngleDegrees;
using ns_cg::Polygon2d;
using ns_cg::SvgStyle;
using ns_cg::Vec2d;
using ns_cg::WriteSvgFile;
using ns_ls2p::CreatePolygonExtractor;
using ns_ls2p::ExtractionMethod;

namespace {

constexpr double kSquareSize = 10.0;
constexpr double kPadding = 1.0;
const Vec2d kCorner(0.0, 0.0);  // the true, sharp corner being analyzed

Polygon2d MakeSquare(double size) {
  return Polygon2d({Vec2d(0.0, 0.0), Vec2d(size, 0.0), Vec2d(size, size),
                     Vec2d(0.0, size)});
}

const char* MethodName(ExtractionMethod method) {
  switch (method) {
    case ExtractionMethod::kMarchingSquares:
      return "MarchingSquares";
    case ExtractionMethod::kMarchingSquaresCornerSharpened:
      return "CornerSharpened";
    case ExtractionMethod::kDualContouring:
      return "DualContouring";
    case ExtractionMethod::kRectilinearThreshold:
      return "RectilinearThreshold";
  }
  return "?";
}

// Distance from `corner` to the outer loop's nearest vertex, and that
// vertex's interior angle.
std::pair<double, double> NearestVertexDistanceAndAngle(
    const std::vector<Vec2d>& outer, const Vec2d& corner) {
  std::size_t best_i = 0;
  double best_d = std::numeric_limits<double>::max();
  for (std::size_t i = 0; i < outer.size(); ++i) {
    const double d = (outer[i] - corner).norm();
    if (d < best_d) {
      best_d = d;
      best_i = i;
    }
  }
  const std::size_t n = outer.size();
  const double angle = InteriorAngleDegrees(outer[(best_i + n - 1) % n],
                                             outer[best_i], outer[(best_i + 1) % n]);
  return {best_d, angle};
}

}  // namespace

int main() {
  const Polygon2d square = MakeSquare(kSquareSize);
  const std::vector<ExtractionMethod> methods = {
      ExtractionMethod::kMarchingSquares,
      ExtractionMethod::kMarchingSquaresCornerSharpened,
      ExtractionMethod::kDualContouring,
      ExtractionMethod::kRectilinearThreshold,
  };

  std::cout << std::fixed << std::setprecision(4);
  std::cout << "True corner: (" << kCorner.x() << ", " << kCorner.y()
            << "), interior angle 90.0000 degrees (exact, by construction)\n\n";

  for (int cells_across : {8, 120}) {
    const std::size_t n = static_cast<std::size_t>(cells_across) + 1;
    const Grid2d<double> field = BuildLevelSet(square, n, n, kPadding);

    std::cout << "cells_across=" << cells_across << " (cell_size=" << field.dx()
              << ")\n";
    std::cout << "  method                | nearest vertex distance | angle\n";
    std::cout << "  ----------------------+-------------------------+---------\n";

    for (ExtractionMethod method : methods) {
      const std::vector<Polygon2d> reconstructed =
          CreatePolygonExtractor(method)->Extract(field);
      const auto [dist, angle] =
          NearestVertexDistanceAndAngle(reconstructed.front().GetOuter(), kCorner);
      std::cout << "  " << std::left << std::setw(22) << MethodName(method)
                 << std::right << " | " << std::setw(23) << dist << " | "
                 << std::setw(8) << angle << "\n";

      if (cells_across == 8) {
        SvgStyle style;
        style.fill = "none";
        style.stroke = "#2563eb";
        style.stroke_width = 0.08;
        const std::string path =
            std::string("square_") + MethodName(method) + ".svg";
        WriteSvgFile(path, reconstructed, style);
      }
    }
    std::cout << "\n";
  }

  SvgStyle sharp_style;
  sharp_style.fill = "none";
  sharp_style.stroke = "#dc2626";
  sharp_style.stroke_width = 0.08;
  WriteSvgFile("square_Original.svg", std::vector<Polygon2d>{square}, sharp_style);

  std::cout << "MarchingSquares stays pinned at a 135-degree angle at both\n"
               "resolutions -- only the cut's size shrinks, never its shape.\n"
               "RectilinearThreshold and CornerSharpened both give (near-)\n"
               "exact 90-degree angles even at the coarse resolution.\n"
               "DualContouring's angle only approaches 90 degrees at the fine\n"
               "resolution (its finite-difference gradient estimate near the\n"
               "corner's singularity is too coarse to place the vertex well\n"
               "at low resolution), but even at the coarse resolution its\n"
               "vertex position is measurably closer to the true corner than\n"
               "MarchingSquares' edge-constrained placement.\n";
  return 0;
}
