#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "common_geometry/levelset.hpp"
#include "common_geometry/math.hpp"
#include "levelset2d_polygon/levelset2d_polygon.hpp"

namespace ns_ls2p {
namespace {

using ns_cg::BuildLevelSet;
using ns_cg::Grid2d;
using ns_cg::Polygon2d;
using ns_cg::Vec2d;

double PolygonArea(const std::vector<Vec2d>& loop) {
  double area = 0.0;
  const std::size_t n = loop.size();
  for (std::size_t i = 0; i < n; ++i) area += ns_cg::Cross2d(loop[i], loop[(i + 1) % n]);
  return 0.5 * area;
}

TEST(MarchingSquaresTest, SquareRoundTrip) {
  Polygon2d square({Vec2d(0.0, 0.0), Vec2d(10.0, 0.0), Vec2d(10.0, 10.0),
                     Vec2d(0.0, 10.0)});
  Grid2d<double> field = BuildLevelSet(square, 41, 41, /*padding=*/2.0);

  const std::vector<Polygon2d> reconstructed = ExtractPolygons(field);

  ASSERT_EQ(reconstructed.size(), 1u);
  EXPECT_TRUE(reconstructed[0].GetHoles().empty());
  EXPECT_NEAR(PolygonArea(reconstructed[0].GetOuter()), 100.0, 1.0);
  EXPECT_GT(PolygonArea(reconstructed[0].GetOuter()), 0.0);  // CCW
}

TEST(MarchingSquaresTest, RingRoundTripHasOneHole) {
  Polygon2d ring(
      {Vec2d(0.0, 0.0), Vec2d(10.0, 0.0), Vec2d(10.0, 10.0), Vec2d(0.0, 10.0)},
      {{Vec2d(3.0, 3.0), Vec2d(3.0, 7.0), Vec2d(7.0, 7.0), Vec2d(7.0, 3.0)}});
  Grid2d<double> field = BuildLevelSet(ring, 61, 61, /*padding=*/2.0);

  const std::vector<Polygon2d> reconstructed = ExtractPolygons(field);

  ASSERT_EQ(reconstructed.size(), 1u);
  EXPECT_NEAR(PolygonArea(reconstructed[0].GetOuter()), 100.0, 1.0);
  EXPECT_GT(PolygonArea(reconstructed[0].GetOuter()), 0.0);  // outer CCW

  ASSERT_EQ(reconstructed[0].GetHoles().size(), 1u);
  const double hole_area = PolygonArea(reconstructed[0].GetHoles()[0]);
  EXPECT_NEAR(std::abs(hole_area), 16.0, 1.0);
  EXPECT_LT(hole_area, 0.0);  // hole CW
}

// Regression test: an off-center, non-square hole sampled at a grid
// resolution that happens to land sample rows exactly on the hole's
// axis-aligned edges (dy=0.2 evenly divides the y=2/y=8 edges here) used to
// produce dozens of degenerate near-zero-area loops instead of one clean
// hole -- caused by SignedDistanceToPolygon returning exact/near-zero
// values, which both destabilized PointInPolygon's ray cast and made
// marching squares interpolate degenerate zero-length segments at those
// grid points.
TEST(MarchingSquaresTest, OffCenterNonSquareHoleGridAlignedWithEdges) {
  Polygon2d ring(
      {Vec2d(0.0, 0.0), Vec2d(20.0, 0.0), Vec2d(20.0, 10.0), Vec2d(0.0, 10.0)},
      {{Vec2d(2.0, 2.0), Vec2d(2.0, 8.0), Vec2d(6.0, 8.0), Vec2d(6.0, 2.0)}});
  Grid2d<double> field = BuildLevelSet(ring, 121, 61, /*padding=*/1.0);

  const std::vector<Polygon2d> reconstructed = ExtractPolygons(field);

  ASSERT_EQ(reconstructed.size(), 1u);
  EXPECT_NEAR(PolygonArea(reconstructed[0].GetOuter()), 200.0, 1.0);
  ASSERT_EQ(reconstructed[0].GetHoles().size(), 1u);
  EXPECT_NEAR(std::abs(PolygonArea(reconstructed[0].GetHoles()[0])), 24.0, 1.0);
}

TEST(MarchingSquaresTest, CircleRoundTripAreaMatchesAnalytic) {
  // Approximate a circle as a many-sided polygon, then round-trip it
  // through a level set.
  constexpr double kRadius = 5.0;
  constexpr int kSides = 64;
  std::vector<Vec2d> circle_outer;
  circle_outer.reserve(kSides);
  for (int i = 0; i < kSides; ++i) {
    const double theta = 2.0 * M_PI * i / kSides;
    circle_outer.push_back(
        Vec2d(kRadius * std::cos(theta), kRadius * std::sin(theta)));
  }
  Polygon2d circle(circle_outer);
  Grid2d<double> field = BuildLevelSet(circle, 81, 81, /*padding=*/1.0);

  const std::vector<Polygon2d> reconstructed = ExtractPolygons(field);

  ASSERT_EQ(reconstructed.size(), 1u);
  EXPECT_TRUE(reconstructed[0].GetHoles().empty());
  const double analytic_area = M_PI * kRadius * kRadius;
  EXPECT_NEAR(PolygonArea(reconstructed[0].GetOuter()), analytic_area, 1.0);
}

}  // namespace
}  // namespace ns_ls2p
