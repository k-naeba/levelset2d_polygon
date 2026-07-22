#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "common_geometry/levelset.hpp"
#include "common_geometry/math.hpp"
#include "levelset2d_polygon/levelset2d_polygon.hpp"

namespace ns_ls2p {
namespace {

using ns_cg::BuildLevelSet;
using ns_cg::Grid2d;
using ns_cg::InteriorAngleDegrees;
using ns_cg::Polygon2d;
using ns_cg::Vec2d;

double PolygonArea(const std::vector<Vec2d>& loop) {
  double area = 0.0;
  const std::size_t n = loop.size();
  for (std::size_t i = 0; i < n; ++i) area += ns_cg::Cross2d(loop[i], loop[(i + 1) % n]);
  return 0.5 * area;
}

Polygon2d MakeSquare(double size) {
  return Polygon2d({Vec2d(0.0, 0.0), Vec2d(size, 0.0), Vec2d(size, size),
                     Vec2d(0.0, size)});
}

Polygon2d MakeRing() {
  return Polygon2d(
      {Vec2d(0.0, 0.0), Vec2d(10.0, 0.0), Vec2d(10.0, 10.0), Vec2d(0.0, 10.0)},
      {{Vec2d(3.0, 3.0), Vec2d(3.0, 7.0), Vec2d(7.0, 7.0), Vec2d(7.0, 3.0)}});
}

// Interior angle of the reconstructed outer loop's vertex nearest `corner`.
double NearestVertexAngle(const std::vector<Vec2d>& outer, const Vec2d& corner) {
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
  return InteriorAngleDegrees(outer[(best_i + n - 1) % n], outer[best_i],
                               outer[(best_i + 1) % n]);
}

TEST(PolygonExtractorTest, AllMethodsReconstructSquareWithReasonableArea) {
  const Polygon2d square = MakeSquare(10.0);
  const Grid2d<double> field = BuildLevelSet(square, 41, 41, 2.0);

  for (auto method : {ExtractionMethod::kMarchingSquares,
                       ExtractionMethod::kMarchingSquaresCornerSharpened,
                       ExtractionMethod::kDualContouring,
                       ExtractionMethod::kRectilinearThreshold}) {
    auto extractor = CreatePolygonExtractor(method);
    const std::vector<Polygon2d> result = extractor->Extract(field);
    ASSERT_EQ(result.size(), 1u) << "method=" << static_cast<int>(method);
    EXPECT_TRUE(result[0].GetHoles().empty());
    // RectilinearThreshold systematically undershoots area (cells only
    // count as inside once their average sample is negative); the others
    // stay close to 100.
    const double expected = method == ExtractionMethod::kRectilinearThreshold ? 90.0 : 100.0;
    EXPECT_NEAR(PolygonArea(result[0].GetOuter()), expected, 10.0)
        << "method=" << static_cast<int>(method);
  }
}

TEST(PolygonExtractorTest, RectilinearThresholdGivesExactRightAngles) {
  const Polygon2d square = MakeSquare(10.0);
  const Grid2d<double> field = BuildLevelSet(square, 9, 9, 1.0);  // coarse

  auto extractor = CreatePolygonExtractor(ExtractionMethod::kRectilinearThreshold);
  const std::vector<Polygon2d> result = extractor->Extract(field);
  ASSERT_EQ(result.size(), 1u);
  // Every axis-aligned corner is an exact right angle, at any resolution --
  // there's no diagonal interpolation to chamfer it.
  EXPECT_DOUBLE_EQ(NearestVertexAngle(result[0].GetOuter(), Vec2d(0.0, 0.0)), 90.0);
}

TEST(PolygonExtractorTest, CornerSharpenedRecoversRightAngleAtCoarseResolution) {
  const Polygon2d square = MakeSquare(10.0);
  const Grid2d<double> field = BuildLevelSet(square, 9, 9, 1.0);  // coarse: cells_across=8

  auto extractor =
      CreatePolygonExtractor(ExtractionMethod::kMarchingSquaresCornerSharpened);
  const std::vector<Polygon2d> result = extractor->Extract(field);
  ASSERT_EQ(result.size(), 1u);
  EXPECT_NEAR(NearestVertexAngle(result[0].GetOuter(), Vec2d(0.0, 0.0)), 90.0, 1e-6);
}

TEST(PolygonExtractorTest, DualContouringConvergesToRightAngleAtFineResolution) {
  const Polygon2d square = MakeSquare(10.0);
  const Grid2d<double> field = BuildLevelSet(square, 121, 121, 1.0);  // fine

  auto extractor = CreatePolygonExtractor(ExtractionMethod::kDualContouring);
  const std::vector<Polygon2d> result = extractor->Extract(field);
  ASSERT_EQ(result.size(), 1u);
  // At fine resolution the finite-difference gradient estimate is accurate
  // enough for the QEF to essentially recover the true corner (unlike
  // marching squares, which stays pinned at a 135-degree chamfer
  // regardless of resolution).
  EXPECT_NEAR(NearestVertexAngle(result[0].GetOuter(), Vec2d(0.0, 0.0)), 90.0, 1.0);
}

TEST(PolygonExtractorTest, DualContouringMovesCornerCloserThanMarchingSquaresEvenAtCoarseResolution) {
  const Polygon2d square = MakeSquare(10.0);
  const Grid2d<double> field = BuildLevelSet(square, 9, 9, 1.0);  // coarse
  const Vec2d corner(0.0, 0.0);

  auto ms = CreatePolygonExtractor(ExtractionMethod::kMarchingSquares)->Extract(field);
  auto dc = CreatePolygonExtractor(ExtractionMethod::kDualContouring)->Extract(field);

  auto nearest_dist = [&](const std::vector<Vec2d>& outer) {
    double best = std::numeric_limits<double>::max();
    for (const auto& v : outer) best = std::min(best, (v - corner).norm());
    return best;
  };
  // The angle itself isn't yet close to 90 at this coarse a resolution
  // (gradient estimation near the corner's singularity is still rough),
  // but the vertex position is already measurably closer to the true
  // corner than marching squares' edge-constrained placement.
  EXPECT_LT(nearest_dist(dc[0].GetOuter()), nearest_dist(ms[0].GetOuter()));
}

TEST(PolygonExtractorTest, CornerSharpenedAndDualContouringDoNotCorruptCircle) {
  constexpr double kRadius = 5.0;
  constexpr int kSides = 64;
  std::vector<Vec2d> circle_outer;
  for (int i = 0; i < kSides; ++i) {
    const double theta = 2.0 * M_PI * i / kSides;
    circle_outer.push_back(
        Vec2d(kRadius * std::cos(theta) + 7.0, kRadius * std::sin(theta) + 7.0));
  }
  const Polygon2d circle(circle_outer);
  const Grid2d<double> field = BuildLevelSet(circle, 81, 81, 1.0);
  const double analytic_area = M_PI * kRadius * kRadius;

  for (auto method : {ExtractionMethod::kMarchingSquaresCornerSharpened,
                       ExtractionMethod::kDualContouring}) {
    auto extractor = CreatePolygonExtractor(method);
    const std::vector<Polygon2d> result = extractor->Extract(field);
    ASSERT_EQ(result.size(), 1u) << "method=" << static_cast<int>(method);
    // A smooth curve has no 135-degree chamfer pairs and no sharp
    // corners, so neither post-process should meaningfully perturb it.
    EXPECT_NEAR(PolygonArea(result[0].GetOuter()), analytic_area, 1.0)
        << "method=" << static_cast<int>(method);
  }
}

TEST(PolygonExtractorTest, AllMethodsHandleRingWithHole) {
  const Polygon2d ring = MakeRing();
  const Grid2d<double> field = BuildLevelSet(ring, 61, 61, 2.0);

  for (auto method : {ExtractionMethod::kMarchingSquares,
                       ExtractionMethod::kMarchingSquaresCornerSharpened,
                       ExtractionMethod::kDualContouring,
                       ExtractionMethod::kRectilinearThreshold}) {
    auto extractor = CreatePolygonExtractor(method);
    const std::vector<Polygon2d> result = extractor->Extract(field);
    ASSERT_EQ(result.size(), 1u) << "method=" << static_cast<int>(method);
    ASSERT_EQ(result[0].GetHoles().size(), 1u) << "method=" << static_cast<int>(method);
    EXPECT_NEAR(PolygonArea(result[0].GetOuter()), 100.0, 5.0)
        << "method=" << static_cast<int>(method);
    EXPECT_NEAR(std::abs(PolygonArea(result[0].GetHoles()[0])), 16.0, 5.0)
        << "method=" << static_cast<int>(method);
  }
}

}  // namespace
}  // namespace ns_ls2p
