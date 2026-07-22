# levelset2d_polygon

A C++17 header-only library that reconstructs a `ns_cg::Polygon2d` from a
2D signed-distance level set (`ns_cg::Grid2d<double>`), via a choice of 4
extraction algorithms (marching squares, marching squares + corner
sharpening, dual contouring, and rectilinear thresholding) selected through
a `PolygonExtractor` factory. Combined with `common_geometry`'s
`BuildLevelSet()`, this closes the round trip: `Polygon2d` -> level set ->
`Polygon2d`.

## Requirements

- CMake >= 3.20
- A C++17 compiler
- [Eigen3](https://eigen.tuxfamily.org/) (e.g. `brew install eigen` on macOS)
- A sibling checkout of [`common_geometry`](../common_geometry) at
  `../common_geometry` relative to this repository

## Building and testing

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build
cd build && ctest --output-on-failure
```

Note: because `common_geometry` is pulled in via `add_subdirectory`, `ctest`
also runs `common_geometry`'s own test suite alongside this project's.

## Usage

```cpp
#include "common_geometry/levelset.hpp"
#include "levelset2d_polygon/levelset2d_polygon.hpp"

ns_cg::Polygon2d ring(/* outer */ {...}, /* holes */ {{...}});
ns_cg::Grid2d<double> field = ns_cg::BuildLevelSet(ring, /*nx=*/81, /*ny=*/81,
                                                    /*padding=*/2.0);

std::vector<ns_cg::Polygon2d> reconstructed = ns_ls2p::ExtractPolygons(field);
// reconstructed[0].GetOuter() / .GetHoles() approximate the original ring.
```

See `examples/roundtrip_demo.cpp` for a runnable version, which also exports
both the original and reconstructed shapes to SVG (via `common_geometry`'s
`svg.hpp`) for visual comparison, plus `levelset_heatmap.svg`: the level
set itself rendered as a diverging-colormap heatmap (blue = inside, red =
outside) with the source polygon outlined on top, via `svg.hpp`'s
`ToSvg(const Grid2d<double>&, const Polygon2d&, ...)` overload.

Its output (regenerate with
`./build/examples/levelset2d_polygon_roundtrip_demo`):

| original | level set heatmap | reconstructed |
| --- | --- | --- |
| <img src="docs/svg/original.svg" width="220"> | <img src="docs/svg/levelset_heatmap.svg" width="220"> | <img src="docs/svg/reconstructed.svg" width="220"> |

## Choosing an extraction method

`ExtractPolygons()` (used above) is always marching squares -- a direct,
dependency-free entry point for the common case. To choose a different
algorithm, go through the `PolygonExtractor` factory instead:

```cpp
#include "levelset2d_polygon/levelset2d_polygon.hpp"

std::unique_ptr<ns_ls2p::PolygonExtractor> extractor =
    ns_ls2p::CreatePolygonExtractor(ns_ls2p::ExtractionMethod::kDualContouring);
std::vector<ns_cg::Polygon2d> reconstructed = extractor->Extract(field);
```

| `ExtractionMethod` | What it does | Tradeoff |
| --- | --- | --- |
| `kMarchingSquares` | Interpolates the zero crossing linearly along each cell edge (see Algorithm below). | Simple and fast, but a right-angle corner is always cut by a single diagonal segment. |
| `kMarchingSquaresCornerSharpened` | Runs marching squares, then a post-process that detects vertex pairs shaped like a diagonal chamfer and collapses them back to the sharp corner their neighboring edges imply. | Cheap to add on top of marching squares; a heuristic, so it only helps shapes that actually have a chamfer-like angle pattern. |
| `kDualContouring` | Places each boundary vertex *inside* its cell (not on an edge), positioned via a least-squares fit (QEF) using the surface normal estimated from the field's gradient at each crossing. | Recovers sharp corners very well *once resolution is fine enough for the gradient estimate to be accurate near the corner*; needs gradient estimation and is more involved than marching squares. |
| `kRectilinearThreshold` | Binarizes the field per cell (no interpolation at all) and traces axis-aligned cell boundaries -- the same technique `rectilinear2d_boolean` uses for its own occupancy grids, reimplemented locally here (Grid2d is already uniform, so no coordinate-compression step is needed, and no dependency on that project's private internals is taken). | Right angles come out exactly right at any resolution; diagonal or curved boundaries come out blocky/staircased instead of smooth, and the binarize threshold shrinks area at coarse resolution. |

See `analysis/corner_chamfer_analysis.cpp` below for a head-to-head comparison.

## Algorithm: marching squares

1. **`MarchCell`**: for each grid cell, classify its 4 corners as inside
   (value < 0) or outside based on the level set, and linearly interpolate
   the zero-crossing point along each edge where the sign changes. This
   yields 0-2 boundary segments per cell (16 standard marching-squares
   cases). Segments are oriented so inside is on the left of the direction
   of travel: outer contours trace CCW, holes CW.
2. **Saddle disambiguation**: the two cases where only diagonally-opposite
   corners are inside are ambiguous from the 4 corner samples alone (are
   the two inside corners connected through the cell center, or separate?).
   This is resolved using the average of the 4 corners as an approximation
   of the center value.
3. **`CollectSegments`** runs `MarchCell` over every cell. Segments from two
   cells sharing a grid edge always interpolate that edge's crossing from
   the same two corner samples in the same order, so shared endpoints are
   bit-identical -- linking needs no epsilon.
4. **`LinkIntoLoops`** links the segment soup into closed loops by
   following each vertex's unique outgoing segment (same technique as
   `rectilinear2d_boolean`'s boundary tracer).
5. Each loop's signed area (shoelace) classifies it as an outer boundary
   (positive/CCW) or a hole (negative/CW); each hole is attached to
   whichever outer polygon contains it (`ns_cg::PointInPolygon`).

## Algorithm: dual contouring

Reuses marching squares' case table verbatim (same 16 cases, same saddle
disambiguation) to decide *which pairs* of edge crossings form each
boundary segment, and its already-validated `LinkIntoLoops`/`SignedArea`
to link segments into loops -- only the *vertex position* assigned to each
segment differs:

1. **`EstimateGradient`**: central-difference gradient at each grid node
   (one-sided at the boundary), an approximation of the outward surface
   normal (the sign convention -- negative inside -- makes the gradient
   point toward increasing/outside values naturally).
2. **`InterpolateNormal`**: for each edge crossing, linearly interpolates
   the two corner gradients using the same parameter `t` as the crossing's
   position, then normalizes -- an approximate normal *at* the crossing.
3. **`SolveQef`**: for a segment's two crossings `(p0,n0)`, `(p1,n1)`,
   solves the 2x2 least-squares system minimizing
   `sum_k (n_k . (x - p_k))^2`, placing the vertex where the two
   crossings' local tangent lines best agree -- which can be *inside* the
   cell, not confined to an edge. Falls back to the midpoint of `p0,p1` if
   the system is near-singular (near-parallel normals, i.e. the boundary is
   locally close to straight -- nothing sharp to resolve there anyway).
4. **`MarchCellDual`** mirrors `MarchCell`'s case switch, computing a
   `DualSegment` (the usual straight-line `Edge2d`, purely for linking, plus
   its `SolveQef`-placed vertex) for each segment instead of just the edge.
5. After linking (on the plain `Edge2d`s), each loop's edges are looked up
   in a `(start,end) -> dual vertex` map to build the final polygon.

This depends on gradient accuracy, which finite differences only estimate
well away from where the *true* gradient is discontinuous -- exactly at a
corner. See "Choosing an extraction method" above for how that plays out in
practice (very good once resolution is fine enough, less impressive at
coarse resolution).

## Algorithm: rectilinear threshold

A cell is "inside" if the average of its 4 corner samples is negative (no
interpolation). From there it's the same occupancy-grid technique
`rectilinear2d_boolean` uses (reimplemented locally, since that project's
`detail::` internals aren't public API and `Grid2d` is already uniform so
the coordinate-compression step isn't needed here anyway): 4-connected
flood-fill labeling, axis-aligned boundary-edge collection tagged with
component id, loop linking, collinear-edge merging, and hole-to-outer
attachment via `PointInPolygon`.

## Algorithm: corner sharpening

A post-process over `ExtractPolygons()`'s output, not a new extraction
technique: `detail::SharpenCorners()` scans a loop's vertices for
consecutive pairs `B,C` whose interior angles are both within a tolerance
(default 5 degrees) of a target angle (default 135, marching squares'
characteristic chamfer angle). Where found, it computes where the
neighboring edges (`A->B` and `C->D`, extended) intersect, and replaces
`B,C` with that single point. Because it only fires on angle patterns that
look like a chamfer, it leaves loops without one (e.g. a smooth curve's
vertices) unchanged -- verified in
`tests/test_polygon_extractor.cpp`'s `CornerSharpenedAndDualContouringDoNotCorruptCircle`.

### A pitfall this hit: exact/near-zero level-set samples

Early versions produced dozens of degenerate, near-zero-area loops instead
of one clean contour whenever the sampling grid happened to place a row (or
column) of nodes essentially *on* the polygon's boundary -- e.g. an
axis-aligned edge at a coordinate the grid spacing divides evenly. Two
compounding problems:

- `PointInPolygon`'s ray-casting test is numerically unstable for points
  essentially on an edge or vertex; neighboring samples a floating-point
  epsilon apart could get opposite inside/outside answers.
- A sample value of exactly (or extremely close to) 0 makes the edge
  interpolation `t = va / (va - vb)` land exactly on that grid point from
  *every* adjacent edge, producing zero-length segments.

`common_geometry::SignedDistanceToPolygon()` now floors near-boundary
magnitudes to a small positive constant (treating them as outside) instead
of returning a value near/at exactly 0, sidestepping both issues. See
`tests/test_marching_squares.cpp`'s
`OffCenterNonSquareHoleGridAlignedWithEdges` for the regression test.

### Limitations

Components that touch another component (or themselves) at a single point
are not supported and may produce an incorrect trace -- the same "bowtie"
caveat documented in `rectilinear2d_boolean`'s `ExtractContours()`.
Correctness beyond the hand-derived case table is validated empirically via
area-convergence tests on known shapes (square, ring, circle approximation)
rather than a formal proof.

**Marching squares cannot reproduce a sharp right-angle corner, at any
resolution.** Each grid cell only ever contributes straight segments, so a
90-degree corner is always replaced by a diagonal cut through whichever
cell it falls in, leaving a vertex that's never actually 90 degrees.
`analysis/corner_chamfer_analysis.cpp` measures this for all 4 extraction
methods, at a coarse (4 cells across) and a fine (120 cells across)
resolution:

| original | MarchingSquares | CornerSharpened | DualContouring | RectilinearThreshold |
| --- | --- | --- | --- | --- |
| <img src="docs/svg/square_Original.svg" width="150"> | <img src="docs/svg/square_MarchingSquares.svg" width="150"> | <img src="docs/svg/square_CornerSharpened.svg" width="150"> | <img src="docs/svg/square_DualContouring.svg" width="150"> | <img src="docs/svg/square_RectilinearThreshold.svg" width="150"> |

(images above are the coarse-resolution reconstructions, deliberately at a
very coarse 4x4 grid so the differences between methods are obvious at a
glance)

Regenerate with `./build/analysis/levelset2d_polygon_corner_chamfer_analysis`;
measured output (a 10x10 square, corner at the origin -- distance and angle
of the reconstructed outer loop's vertex nearest that corner):

```
cells_across=4 (cell_size=3.0000)
  method                | nearest vertex distance | angle
  ----------------------+-------------------------+---------
  MarchingSquares        |                  2.0000 | 135.0000
  CornerSharpened        |                  0.0000 |  90.0000
  DualContouring         |                  0.5751 | 118.6357
  RectilinearThreshold   |                  2.2361 |  90.0000

cells_across=120 (cell_size=0.1000)
  method                | nearest vertex distance | angle
  ----------------------+-------------------------+---------
  MarchingSquares        |                  0.1000 | 135.0000
  CornerSharpened        |                  0.0000 |  90.0000
  DualContouring         |                  0.0000 |  90.0000
  RectilinearThreshold   |                  0.0000 |  90.0000
```

- **MarchingSquares** stays pinned at exactly 135 degrees at both
  resolutions -- refining the grid shrinks the cut, never its shape.
- **CornerSharpened** and **RectilinearThreshold** both recover the exact
  right angle even at the coarse resolution (for different reasons: one by
  detecting and undoing the chamfer, the other by never creating diagonal
  edges in the first place).
- **DualContouring** only closes in on exactly 90 degrees at the fine
  resolution -- at the coarse resolution its finite-difference gradient
  estimate near the corner's singularity isn't accurate enough to place the
  vertex perfectly (118.6 degrees here; how far off varies with exactly
  where the corner falls within its cell, so this number isn't necessarily
  monotonic across arbitrary resolutions -- see the git history of this
  file for a run that measured 149 degrees at a different coarse
  resolution). What *is* consistent: its vertex *position* is measurably
  closer to the true corner than marching squares' at every resolution
  tested (0.575 vs. 2.0 here), and it's the only method whose accuracy
  keeps improving as resolution increases rather than staying fixed.

## Directory layout

```
levelset2d_polygon/
├── include/levelset2d_polygon/  Public headers (header-only library)
├── examples/                    Runnable demo
├── analysis/                    Programs characterizing algorithmic behavior
│                                 (not usage examples), e.g. the corner-chamfer
│                                 measurement above
└── tests/                       GoogleTest unit tests
```
