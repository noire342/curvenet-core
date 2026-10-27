# Core architecture and data contracts

`Solver` owns binding, projection, cut-face construction, differential operators,
transport and sparse factorization. `RigSession` validates and caches a complete
binding and maps animator offsets to posed controls. Both are independent of any
host application. Each session owns its solver; calls that mutate or solve through
the same session must be serialized by the host.

## Mesh and network

`Mesh` has finite local-space points and consistently wound manifold polygon
faces. `Net` is a point pool and cubic spans indexed `[knot, handle, handle,
knot]`. Profiles that meet must share knot indices; ambiguous crossings and
overlapping projected curves are rejected. Endpoints may have arbitrary valence.

The projection mesh defines the bind surface. On evaluation, the upstream rest
pose may change positions but must retain its vertex count, polygon connectivity
and index order. Projection, upstream and curve points are in the same local
space. The output retains the input vertex topology.

`RigSession` requires one parent entry per control: `-1` for a knot, otherwise
the index of that handle's knot. Handles belong to their span's corresponding
endpoint. Knots expose translation and a 3x3 linear transform; handles expose
translation only. Inputs must be finite.

For warped base controls `b`, transported frames `R`, knot transform `T` and
offsets `d`, the knot is `b_i + R_i d_i`; its handle is
`b_i + R_i [T_i(R_i^T(b_h-b_i) + d_h) + d_i]`.

## Binding lifecycle

`ensureBinding` compares projection coordinates and topology, all network points,
spans, parent indices, options and explicit revision. A change reconstructs the
binding; invalid changes are rejected before they can index an old descriptor.
Neutral channels pass upstream positions through exactly. Moving a knot changes
the associated cubic handles and surface constraints; the factorization is reused
while the binding remains unchanged.

## Numerical approach and limits

The implementation follows the profile-curve paper's broad construction: sample
curves, project to the surface, construct local cut faces, use polygon operators
and solve sparse systems. Eigen provides the linear algebra and SimplicialLDLT.
The cutter uses bounded unfolded corridors, a short adjacent-face edge fallback,
and a shared-vertex fallback at singularities. It does not establish equivalence
to Pixar's internal implementation. Sparse residuals measure numerical solution
quality, not anatomical fidelity or the visual appeal of an articulated pose.

The core does not implement a host dependency graph, user interface, curve layout
generation, skeleton binding, collision handling, self-contact, or mesh export.
Consumers must arrange those separately and validate their own coordinate and
lifecycle contracts.
