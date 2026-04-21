# Repository Guidelines

## Project Structure & Module Organization

This repository is a cross-platform C++20 CLI that converts SplatKing / ARKit LiDAR captures into COLMAP models.

- `src/`: main application code.
  - `main.cpp`: CLI entrypoint.
  - `sk2cm.hpp`: shared types and public API.
  - `sk2cm.cpp`: conversion pipeline, I/O, validation, and alignment logic.
- `third_party/nlohmann/json.hpp`: vendored JSON dependency. Do not modify unless intentionally upgrading it.
- `build/`: local CMake build output. Ignored by Git.

There is currently no dedicated `tests/` directory.

## Build, Test, and Development Commands

- `cmake -S . -B build`
  Configures the project.
- `cmake --build build -j4`
  Builds the `sk2colmap` executable.
- `./build/sk2colmap --help`
  Prints CLI usage and confirms the binary runs.
- `./build/sk2colmap --input ./LidarSeries --output ./output --validate-exported-model`
  Runs a full local export and checks generated COLMAP text/bin consistency.

When validating behavior changes, prefer a real dataset and inspect:
- `output/model_validation.txt`
- `output/validation_summary.txt`
- `output/reprojection_summary.txt`

## Coding Style & Naming Conventions

- Use C++20.
- Use 4-space indentation and keep brace style consistent with existing files.
- Prefer clear, descriptive names over abbreviations.
- Types use `PascalCase` (`LoadedDataset`), functions use `PascalCase` for public API and local helpers matching existing code, variables use `snake_case` or descriptive lower-case names consistent with surrounding code.
- Keep new code in ASCII unless a file already requires otherwise.

## Testing Guidelines

There is no automated unit test suite yet. Minimum validation for changes:

1. Rebuild with CMake.
2. Run `./build/sk2colmap --help`.
3. Run at least one real conversion with `--validate-exported-model`.
4. If touching geometry, depth mapping, or alignment, compare output counts and summaries against a known-good dataset.

## Commit & Pull Request Guidelines

Commits should use short imperative subjects, e.g. `Port SplatKing converter to cross-platform C++`.

For pull requests:
- explain the behavior change,
- list validation commands you ran,
- call out any output-count or numeric differences,
- include sample paths or repro steps for dataset-dependent issues.

## Data & Repository Hygiene

- Do not commit datasets, generated `output/` folders, or CMake build artifacts.
- Keep large binary assets out of the repository unless they are essential release artifacts.
