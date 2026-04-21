# SK2CM

Native C++ converter for exporting SplatKing / ARKit LiDAR captures to a COLMAP model. The binary targets Linux, macOS, and Windows via CMake.

## Prerequisites

- CMake 3.20+
- C++20 compiler
  - Linux: GCC 11+ or Clang 14+
  - macOS: recent Apple Clang
  - Windows: Visual Studio 2022 or recent `clang-cl`

## Build

From the repository root:

```bash
cmake -S . -B build
cmake --build build --config Release
```

Produced executables:

- Linux / macOS: `build/sk2colmap`
- Windows multi-config: `build/Release/sk2colmap.exe`

## Usage

Example:

```bash
./build/sk2colmap --input ./LidarSeries --output ./output_final --validate-exported-model
```

Useful options:

- `--align-to-colmap-model <Path>` aligns the export to the coordinate frame of a reference COLMAP model
- `--snap-shared-cameras-to-reference` forces shared image poses to match the reference model exactly
- `--no-copy-images` skips copying JPEGs into the output directory
- `--validate-exported-model` rereads the exported `.txt` and `.bin` files and compares their counts
- `--validation-images <N>` changes how many images are used for depth validation
- `--validation-stride <N>` changes LiDAR point subsampling during validation
- `--depth-bucket-size <N>` controls spatial subsampling during LiDAR projection
- `--max-observations-per-track <N>` limits how many observations are kept for each 3D point
- `--min-capture-gap <N>` enforces a minimum frame gap between observations in the same track
- `--min-baseline-m <meters>` and `--min-angle-deg <degrees>` require a minimum geometric baseline
- `--base-depth-tolerance-m <meters>` and `--relative-depth-tolerance <ratio>` tune depth filtering

## Outputs

The converter writes:

- `output/sparse/0/cameras.bin`
- `output/sparse/0/cameras.txt`
- `output/sparse/0/images.bin`
- `output/sparse/0/images.txt`
- `output/sparse/0/points3D.bin`
- `output/sparse/0/points3D.txt`
- `output/images/`
- `output/lidar_dense_world_points.ply`
- `output/lidar_sparse_world_points.ply`
- `output/model_validation.txt` quand `--validate-exported-model` est active
- `output/alignment_summary.txt` quand `--align-to-colmap-model` est utilise
- `output/camera_pose_summary.txt`
- `output/validation_summary.txt`
- `output/reprojection_summary.txt`

## Current Behavior

- Poses and intrinsics are exported in both COLMAP text and binary formats.
- Portrait JPEGs stored as `1440x1920` are automatically remapped from the original landscape intrinsics and depth maps.
- A separate `camera_id` is created per image because intrinsics vary slightly during capture.
- `lidar_dense_world_points.ply` is the raw session-level SplatKing point cloud.
- `lidar_sparse_world_points.ply` and `points3D` are built from real LiDAR points, keeping only image/depth-consistent observations.
- `--snap-shared-cameras-to-reference` prioritizes visual alignment of shared cameras with a reference COLMAP model, with a possible tradeoff in LiDAR reprojection consistency.
- `camera_pose_summary.txt` lets you numerically verify that exported camera centers match the source ARKit poses.
- Track refinement removes overly redundant observations and requires a minimum geometric baseline between views.

## Notes

- The portable implementation lives in `src/` with a root-level CMake build.
- Recent validation on `LidarSeries_20260420_105850` matched the C# version on camera, image, and 3D point counts, with one additional 2D observation retained on the C++ side due to double-precision arithmetic and a round-to-even fix in depth mapping.
