// Demonstrates the full round trip: Polygon2d -> level set (Grid2d<double>)
// -> Polygon2d (via marching squares). Uses a ring-shaped polygon (a square
// with a square hole) so both outer-contour and hole reconstruction are
// exercised. Exports both the original and reconstructed shapes to SVG for
// visual comparison.

#include <iostream>

#include "common_geometry/levelset.hpp"
#include "common_geometry/svg.hpp"
#include "levelset2d_polygon/levelset2d_polygon.hpp"

using ns_cg::BuildLevelSet;
using ns_cg::Grid2d;
using ns_cg::Polygon2d;
using ns_cg::SvgStyle;
using ns_cg::Vec2d;
using ns_cg::WriteSvgFile;

namespace {

double PolygonArea(const std::vector<Vec2d>& loop) {
  double area = 0.0;
  const std::size_t n = loop.size();
  for (std::size_t i = 0; i < n; ++i)
    area += ns_cg::Cross2d(loop[i], loop[(i + 1) % n]);
  return 0.5 * area;
}

}  // namespace

int main() {
  const Polygon2d original(
      {Vec2d(0.0, 0.0), Vec2d(10.0, 0.0), Vec2d(10.0, 10.0), Vec2d(0.0, 10.0)},
      {{Vec2d(3.0, 3.0), Vec2d(3.0, 7.0), Vec2d(7.0, 7.0), Vec2d(7.0, 3.0)}});

  const Grid2d<double> field = BuildLevelSet(original, 81, 81, /*padding=*/2.0);
  const std::vector<Polygon2d> reconstructed = ns_ls2p::ExtractPolygons(field);

  std::cout << "original: outer=" << original.GetOuter().size()
            << " vertices, area=" << PolygonArea(original.GetOuter())
            << ", holes=" << original.GetHoles().size() << "\n";
  for (const auto& poly : reconstructed) {
    std::cout << "reconstructed: outer=" << poly.GetOuter().size()
              << " vertices, area=" << PolygonArea(poly.GetOuter())
              << ", holes=" << poly.GetHoles().size();
    for (const auto& hole : poly.GetHoles())
      std::cout << " (hole area=" << PolygonArea(hole) << ")";
    std::cout << "\n";
  }

  // GitHub's README renderer does not reliably honor an <img> tag's
  // width/height override for locally-referenced SVGs, so the desired
  // display size is baked into the SVG itself here (via target_width_px)
  // rather than relied on at the markdown layer -- see the README's
  // image tables, which reference these files with no size attributes.
  SvgStyle style;
  style.fill = "#22c55e";
  style.stroke = "black";
  style.stroke_width = 0.05;
  style.target_width_px = 200.0;
  WriteSvgFile("original.svg", std::vector<Polygon2d>{original}, style);
  WriteSvgFile("reconstructed.svg", reconstructed, style);

  SvgStyle heatmap_style;
  heatmap_style.stroke = "black";
  heatmap_style.stroke_width = 0.1;
  heatmap_style.target_width_px = 200.0;
  WriteSvgFile("levelset_heatmap.svg", field, original, heatmap_style);

  std::cout << "wrote original.svg, reconstructed.svg, levelset_heatmap.svg\n";
  return 0;
}
