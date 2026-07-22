#pragma once

#include <cstddef>
#include <vector>

#include "common_geometry/grid.hpp"
#include "common_geometry/polygon.hpp"
#include "common_geometry/rectilinear_trace.hpp"
#include "common_geometry/types.hpp"

// Internal implementation of the "rectilinear threshold" extraction method:
// binarize the level set onto the grid's own cells (no interpolation) and
// trace axis-aligned boundaries via common_geometry's shared
// TraceRectilinearBoundary() -- the same technique rectilinear2d_boolean
// uses for its own (unrelated, BBox2d-set-derived) occupancy grids. Not
// part of the public API.
namespace ns_ls2p::detail {

// Cell (col,row) [corners at grid nodes (col,row)..(col+1,row+1)] is
// "inside" if the average of its 4 corner samples is negative.
inline bool CellInside(const ns_cg::Grid2d<double>& field, std::size_t col,
                        std::size_t row) {
  const double avg = (field.at(col, row) + field.at(col + 1, row) +
                       field.at(col + 1, row + 1) + field.at(col, row + 1)) *
                      0.25;
  return avg < 0.0;
}

inline std::vector<ns_cg::Polygon2d> ExtractPolygonsRectilinearThreshold(
    const ns_cg::Grid2d<double>& field) {
  if (field.nx() < 2 || field.ny() < 2) return {};

  const std::size_t num_cols = field.nx() - 1;
  const std::size_t num_rows = field.ny() - 1;

  std::vector<char> occupancy(num_cols * num_rows);
  for (std::size_t row = 0; row < num_rows; ++row)
    for (std::size_t col = 0; col < num_cols; ++col)
      occupancy[row * num_cols + col] = CellInside(field, col, row) ? 1 : 0;

  return ns_cg::TraceRectilinearBoundary(
      occupancy, num_cols, num_rows,
      [&](std::size_t col, std::size_t row) { return field.position(col, row); });
}

}  // namespace ns_ls2p::detail
