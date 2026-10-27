# CurveNet Core

Experimental C++17 surface deformation driven by a network of cubic Bézier
profile curves. This repository contains the host-independent numerical core and
its tests, plus an optional Maya 2023 C++ deformer adapter and authoring helper.
The core alone does not require Maya or Blender.

The Maya 2023 adapter now embeds its authoring helper in a single `.mll`.
Load it, select one polygon mesh and one or more NURBS curves, then choose
**CurveNet > Bind selected mesh and curves** (or `maya.cmds.cnBindSelection()`).
See the [Maya workflow and limitations](docs/maya2023.md).

The implementation is inspired by de Goes, Sheffler and Fleischer,
[Character Articulation through Profile Curves (2022)](https://research.pixar.com/docs/2022.SiggraphPapers.GSF.pdf).
It is an independent research implementation, not Pixar software. Equivalence to
Pixar's production system or production-ready deformation quality is not established.

## What is here

- `cn::Solver`: polygon meshes, Bézier graph, projection, cut-aware assembly,
  surface frame transport, sparse factorization and repeated solves.
- `cn::RigSession`: cache invalidation, topology and parent validation, and
  direct knot/tangent channel mapping over an upstream pose.
- `tests/`: synthetic geometry tests covering invariants and failure cases.
- `examples/basic.cpp`: a small executable showing binding and deformation.
- `docs/architecture.md`: spaces, graph format, cache and limitations.
- `maya/`: optional Maya 2023 `MPxDeformerNode`, Python setup and smoke test.
- [Maya build and use](docs/maya2023.md): create a demo or bind multiple curves.

## Build and test

Requirements: CMake 3.24+, a C++17 compiler, and Eigen 3.4.0. CMake searches for
an installed Eigen first, then downloads Eigen 3.4.0 with a SHA-256 check. For an
offline build, provide the directory containing `Eigen/Core` with `EIGEN_ROOT`.

```powershell
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

```powershell
cmake -S . -B build -DEIGEN_ROOT=C:/path/to/eigen-3.4.0
```

The target `CurveNet::core` exposes Eigen's C++ types. The core library contains no
Autodesk or Blender SDK code; `maya/` provides a separate optional adapter.
Keep a binding and use it for multiple solves while projection geometry, graph and
options remain unchanged. The public API is experimental and may change.

## Status

Synthetic host tests cover neutral poses, rigid transformations, cuts, junctions,
upstream changes, cache invalidation, direct controls and malformed data. Those
tests do not establish realistic anatomy, animation quality, self-contact, Maya
runtime compatibility, or performance on arbitrary production assets.

The experimental cutter rejects ambiguous overlaps and unsupported projection
jumps. See [architecture and limitations](docs/architecture.md).

## License

The project is licensed under [GPL-3.0-or-later](LICENSE). Eigen is a separate
dependency and its [MPL-2.0 license](docs/EIGEN-LICENSE.txt) applies to Eigen.
The public repository does not contain the reference paper, Pixar assets, or user
models. All trademarks belong to their respective owners.
