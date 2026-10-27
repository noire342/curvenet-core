"""Maya 2023 authoring adapter; all deformation runs in the C++ plug-in.

No NumPy, callbacks or script nodes are required. Channels are mesh-local offsets
in the transported projection frame. Tangent handles expose translation only.
"""

import json
import math
import itertools
from pathlib import Path


def validate_bundle(data):
    """Validate the mesh-local interchange format without loading Maya."""
    if data.get("schema") != "curvenet.bundle/1" or data.get("space") != "mesh-local":
        raise ValueError("Expected curvenet.bundle/1 in mesh-local coordinates")
    projection, net, frames = data["projection"], data["network"], data["frames"]

    def vectors(values, count=None):
        if count is not None and len(values) != count:
            raise ValueError("Array length mismatch")
        if not values or any(
            len(p) != 3
            or not all(isinstance(v, (int, float)) and math.isfinite(v) for v in p)
            for p in values
        ):
            raise ValueError("Expected finite 3D vectors")

    vectors(projection["points"])
    vectors(net["points"])
    count = len(net["points"])
    parents = net["parents"]
    if len(parents) != count:
        raise ValueError("Parent count mismatch")
    for i, parent in enumerate(parents):
        if (
            type(parent) is not int
            or parent < -1
            or parent >= count
            or (parent >= 0 and (parent == i or parents[parent] != -1))
        ):
            raise ValueError("Invalid tangent parent")
    if not net["splines"]:
        raise ValueError("Network has no splines")
    for spline in net["splines"]:
        if len(spline) != 4 or any(
            type(i) is not int or i < 0 or i >= count for i in spline
        ):
            raise ValueError("Invalid spline indices")
        a, h1, h2, b = spline
        if parents[a] != -1 or parents[b] != -1 or parents[h1] != a or parents[h2] != b:
            raise ValueError("Spline handles must belong to their endpoint knots")
    for face in projection["faces"]:
        if (
            len(face) < 3
            or len(set(face)) != len(face)
            or any(
                type(i) is not int or i < 0 or i >= len(projection["points"])
                for i in face
            )
        ):
            raise ValueError("Invalid mesh face")
    vectors(frames["positions"], count)
    if len(frames["rotations"]) != count:
        raise ValueError("Frame count mismatch")
    for rotation in frames["rotations"]:
        vectors(rotation, 3)
        for i in range(3):
            for j in range(3):
                expected = 1.0 if i == j else 0.0
                if (
                    abs(
                        sum(rotation[k][i] * rotation[k][j] for k in range(3))
                        - expected
                    )
                    > 1e-5
                ):
                    raise ValueError("Control frame must be orthonormal")
    return data


def read_bundle(path):
    """Read JSON data; never execute scripts from an interchange file."""
    return validate_bundle(json.loads(Path(path).read_text(encoding="utf-8")))


def load_plugin(path=None):
    """Load a Maya 2023 binary, or reuse the registered cnProfileMover type."""
    from maya import cmds

    if int(cmds.about(apiVersion=True)) // 10000 != 2023:
        raise RuntimeError("This adapter requires Maya 2023")
    for plugin in cmds.pluginInfo(query=True, listPlugins=True) or []:
        if "cnProfileMover" in (
            cmds.pluginInfo(plugin, query=True, dependNode=True) or []
        ):
            if path is not None:
                loaded = Path(cmds.pluginInfo(plugin, query=True, path=True)).resolve()
                if loaded != Path(path).resolve():
                    raise RuntimeError(
                        "Another CurveNet binary is already loaded: %s. Restart Maya before loading %s."
                        % (loaded, path)
                    )
            return plugin
    plugin = cmds.loadPlugin(str(path) if path is not None else "curvenet", quiet=True)[
        0
    ]
    if "cnProfileMover" not in (
        cmds.pluginInfo(plugin, query=True, dependNode=True) or []
    ):
        raise RuntimeError("The loaded plug-in did not register cnProfileMover")
    return plugin


def _matrix(rotation, translation):
    # Eigen/JSON column-vector convention -> Maya row-vector convention.
    return [
        rotation[j][i]
        if i < 3 and j < 3
        else translation[j]
        if i == 3 and j < 3
        else 1.0
        if i == j
        else 0.0
        for i in range(4)
        for j in range(4)
    ]


def create_rig(data, mesh=None, make_controls=True):
    """Create a test mesh, or attach after the existing history of `mesh`.

    For an existing rig/cache, topology and vertex order must match the bundle.
    The projection snapshot and current upstream shape may have different poses.
    The source_world_matrix field is metadata, never an implicit axis conversion.
    """
    from maya import cmds
    from maya.api import OpenMaya as om

    data = validate_bundle(data)
    load_plugin()
    projection, net = data["projection"], data["network"]
    counts = [len(face) for face in projection["faces"]]
    indices = [i for face in projection["faces"] for i in face]

    def mesh_data():
        owner = om.MFnMeshData().create()
        om.MFnMesh().create(
            [om.MPoint(*p) for p in projection["points"]], counts, indices, parent=owner
        )
        return owner

    def obj(name):
        selection = om.MSelectionList()
        selection.add(name)
        return selection.getDependNode(0)

    def validate_mesh(mesh):
        shapes = (
            cmds.listRelatives(mesh, shapes=True, noIntermediate=True, fullPath=True)
            or []
        )
        if len(shapes) != 1 or cmds.nodeType(shapes[0]) != "mesh":
            raise ValueError("Choose a transform with exactly one polygon shape")
        incoming = om.MFnMesh(obj(shapes[0]))
        actual_counts, actual_indices = incoming.getVertices()
        if (
            incoming.numVertices != len(projection["points"])
            or list(actual_counts) != counts
            or list(actual_indices) != indices
        ):
            raise ValueError(
                "Topology or vertex order differs from the exported projection"
            )

    if mesh is not None:
        validate_mesh(mesh)  # Never undo a previous user action on invalid input.
    mutation_started = False
    cmds.undoInfo(openChunk=True, chunkName="Import CurveNet")
    try:
        if mesh is None:
            mesh = cmds.createNode("transform", name="CurveNetBody")
            mutation_started = True
            shape = cmds.createNode("mesh", name="CurveNetBodyShape", parent=mesh)
            om.MFnDependencyNode(obj(shape)).findPlug("inMesh", False).setMObject(
                mesh_data()
            )
            cmds.sets(shape, edit=True, forceElement="initialShadingGroup")
        node = cmds.deformer(
            mesh, type="cnProfileMover", name="CurveNetDeformer", after=True
        )[0]
        mutation_started = True
        cmds.setAttr(node + ".envelope", 0)
        fn = om.MFnDependencyNode(obj(node))
        fn.findPlug("projectionMesh", False).setMObject(mesh_data())
        fn.findPlug("restControlPoints", False).setMObject(
            om.MFnPointArrayData().create([om.MPoint(*p) for p in net["points"]])
        )
        fn.findPlug("splineIndices", False).setMObject(
            om.MFnIntArrayData().create([i for sp in net["splines"] for i in sp])
        )
        fn.findPlug("controlParents", False).setMObject(
            om.MFnIntArrayData().create(net["parents"])
        )
        controls, frames = {}, data["frames"]
        size = (
            max(
                max(p[k] for p in projection["points"])
                - min(p[k] for p in projection["points"])
                for k in range(3)
            )
            * 0.025
        )
        size = max(size, 0.01)
        group = None
        if make_controls:
            group = cmds.createNode("transform", name="CurveNetControls", parent=mesh)
            for i, parent in enumerate(net["parents"]):
                if parent >= 0:
                    continue
                anchor = cmds.createNode(
                    "transform", name="CN_Frame_%04d" % i, parent=group
                )
                cmds.xform(
                    anchor,
                    objectSpace=True,
                    matrix=_matrix(frames["rotations"][i], frames["positions"][i]),
                )
                control = cmds.createNode(
                    "transform", name="CN_Knot_%04d" % i, parent=anchor
                )
                locator = cmds.createNode("locator", parent=control)
                cmds.setAttr(locator + ".localScale", size, size, size)
                cmds.setAttr(locator + ".overrideEnabled", True)
                cmds.setAttr(locator + ".overrideColor", 9)
                cmds.connectAttr(
                    control + ".matrix", "%s.controlOffsets[%d]" % (node, i)
                )
                controls[i] = control
            for i, parent in enumerate(net["parents"]):
                if parent < 0:
                    continue
                delta = [
                    frames["positions"][i][j] - frames["positions"][parent][j]
                    for j in range(3)
                ]
                rotation = frames["rotations"][parent]
                local = [
                    sum(rotation[k][j] * delta[k] for k in range(3)) for j in range(3)
                ]
                anchor = cmds.createNode(
                    "transform", name="CN_HandleFrame_%04d" % i, parent=controls[parent]
                )
                cmds.setAttr(anchor + ".translate", *local)
                control = cmds.createNode(
                    "transform", name="CN_Handle_%04d" % i, parent=anchor
                )
                locator = cmds.createNode("locator", parent=control)
                cmds.setAttr(
                    locator + ".localScale", size * 0.6, size * 0.6, size * 0.6
                )
                cmds.setAttr(locator + ".overrideEnabled", True)
                cmds.setAttr(locator + ".overrideColor", 14)
                for channel in ("rx", "ry", "rz", "sx", "sy", "sz"):
                    cmds.setAttr(control + "." + channel, lock=True, keyable=False)
                cmds.connectAttr(
                    control + ".matrix", "%s.controlOffsets[%d]" % (node, i)
                )
                controls[i] = control
        cmds.setAttr(node + ".envelope", 1)
        cmds.select(mesh)
        result = {"mesh": mesh, "deformer": node, "controls": controls, "group": group}
    except Exception:
        cmds.undoInfo(closeChunk=True)
        if mutation_started:
            cmds.undo()
        raise
    else:
        cmds.undoInfo(closeChunk=True)
    return result


def import_bundle(path, mesh=None, make_controls=True):
    """Create a complete Maya deformer and controls from a JSON bundle."""
    return create_rig(read_bundle(path), mesh=mesh, make_controls=make_controls)


def graph_from_segments(segments, weld_tolerance=1e-4):
    """Convert cubic spans into a shared-knot graph in mesh-local space."""
    if not math.isfinite(weld_tolerance) or weld_tolerance <= 0:
        raise ValueError("Weld tolerance must be positive and finite")
    points, parents, splines, buckets = [], [], [], {}

    def key(p):
        return tuple(math.floor(x / weld_tolerance) for x in p)

    def knot(p):
        cell = key(p)
        for offset in itertools.product((-1, 0, 1), repeat=3):
            adjacent = tuple(cell[i] + offset[i] for i in range(3))
            for index in buckets.get(adjacent, ()):
                if math.dist(points[index], p) <= weld_tolerance:
                    return index
        index = len(points)
        points.append(list(p))
        parents.append(-1)
        buckets.setdefault(cell, []).append(index)
        return index

    for cp in segments:
        if len(cp) != 4 or any(
            len(p) != 3 or not all(math.isfinite(x) for x in p) for p in cp
        ):
            raise ValueError("Each span requires four finite 3D points")
        a, b = knot(cp[0]), knot(cp[3])
        if a == b:
            raise ValueError(
                "A curve span has coincident endpoints; split or edit the curve"
            )
        h = len(points)
        points.extend([list(cp[1]), list(cp[2])])
        parents.extend([a, b])
        splines.append([a, h, h + 1, b])
    if not splines:
        raise ValueError("Select at least one nonempty curve")
    return {"schema": 1, "points": points, "parents": parents, "splines": splines}


def _bezier_from_samples(samples):
    """Recover a cubic polynomial from positions at 0, 1/3, 2/3 and 1."""
    p0, q1, q2, p3 = samples
    u = [27 * q1[i] - 8 * p0[i] - p3[i] for i in range(3)]
    v = [27 * q2[i] - p0[i] - 8 * p3[i] for i in range(3)]
    return [
        p0,
        [(2 * u[i] - v[i]) / 18 for i in range(3)],
        [(2 * v[i] - u[i]) / 18 for i in range(3)],
        p3,
    ]


def bind_curves(mesh, curves, weld_tolerance=1e-4):
    """Bind one mesh to multiple Maya NURBS curves; create direct controls.

    Curve geometry is sampled at bind time; edit the generated controls after
    binding. Cubic/linear polynomial spans are converted exactly. Crossings
    must be shared endpoints or existing span knots. Rational spans rejected.
    """
    from maya import cmds
    from maya.api import OpenMaya as om

    def shape_of(name, kind):
        if cmds.objectType(name, isAType=kind):
            return (cmds.ls(name, long=True) or [name])[0]
        shapes = (
            cmds.listRelatives(
                name, shapes=True, noIntermediate=True, fullPath=True, type=kind
            )
            or []
        )
        if len(shapes) != 1:
            raise ValueError("Expected exactly one %s shape on %s" % (kind, name))
        return shapes[0]

    if isinstance(curves, str):
        curves = [curves]
    if not curves:
        raise ValueError("Supply at least one NURBS curve")
    mesh_shape = shape_of(mesh, "mesh")
    mesh_transform = (
        cmds.listRelatives(mesh_shape, parent=True, fullPath=True) or [None]
    )[0]
    sel = om.MSelectionList()
    sel.add(mesh_shape)
    mesh_path = sel.getDagPath(0)
    mesh_fn = om.MFnMesh(mesh_path)
    inverse = mesh_path.inclusiveMatrixInverse()
    mesh_points = [[p.x, p.y, p.z] for p in mesh_fn.getPoints(om.MSpace.kObject)]
    counts, ids = mesh_fn.getVertices()
    faces, pos = [], 0
    for count in counts:
        faces.append(list(ids[pos : pos + count]))
        pos += count
    scale = max(
        max(p[i] for p in mesh_points) - min(p[i] for p in mesh_points)
        for i in range(3)
    )
    tolerance = max(scale * 1e-8, 1e-8)
    spans = []
    for curve in curves:
        curve_shape = shape_of(curve, "nurbsCurve")
        sel = om.MSelectionList()
        sel.add(curve_shape)
        fn = om.MFnNurbsCurve(sel.getDagPath(0))
        if fn.degree not in (1, 3):
            raise ValueError("Curve %s must have degree 1 or 3" % curve)
        lo, hi = sorted(fn.knotDomain)
        knots = sorted(set([lo, hi] + [float(k) for k in fn.knots() if lo < k < hi]))

        def sample(t):
            point = fn.getPointAtParam(t, om.MSpace.kWorld) * inverse
            return [point.x, point.y, point.z]

        for start, end in zip(knots, knots[1:]):
            if end - start <= 1e-12:
                continue
            cp = _bezier_from_samples(
                [sample(start + (end - start) * i / 3) for i in range(4)]
            )
            for t in (0.2, 0.5, 0.8):
                expected = sample(start + (end - start) * t)
                s = 1 - t
                actual = [
                    s**3 * cp[0][i]
                    + 3 * s * s * t * cp[1][i]
                    + 3 * s * t * t * cp[2][i]
                    + t**3 * cp[3][i]
                    for i in range(3)
                ]
                if math.dist(actual, expected) > tolerance:
                    raise ValueError(
                        "Curve %s contains rational/non-polynomial spans; rebuild it as a non-rational cubic"
                        % curve
                    )
            spans.append(cp)
    graph = graph_from_segments(spans, weld_tolerance)
    identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    bundle = {
        "schema": "curvenet.bundle/1",
        "space": "mesh-local",
        "projection": {"points": mesh_points, "faces": faces},
        "network": graph,
        "frames": {
            "positions": graph["points"],
            "rotations": [identity for _ in graph["points"]],
        },
    }
    return create_rig(bundle, mesh=mesh_transform)


def bind_selection(weld_tolerance=1e-4):
    """Select one mesh and one or more curves, then call this function."""
    from maya import cmds

    meshes, curves = [], []
    for item in cmds.ls(selection=True, long=True) or []:
        shape_types = [cmds.nodeType(item)] + [
            cmds.nodeType(s)
            for s in cmds.listRelatives(
                item, shapes=True, noIntermediate=True, fullPath=True
            )
            or []
        ]
        if "mesh" in shape_types:
            meshes.append(item)
        elif "nurbsCurve" in shape_types:
            curves.append(item)
    if len(meshes) != 1 or not curves:
        raise ValueError("Select exactly one polygon mesh and one or more NURBS curves")
    return bind_curves(meshes[0], curves, weld_tolerance)


def demo_bundle():
    """An original synthetic quad patch and three joined cubic profiles."""
    n = 10
    points = [[x - 5.0, float(y), 0.0] for y in range(n + 1) for x in range(n + 1)]
    faces = [
        [
            y * (n + 1) + x,
            y * (n + 1) + x + 1,
            (y + 1) * (n + 1) + x + 1,
            (y + 1) * (n + 1) + x,
        ]
        for y in range(n)
        for x in range(n)
    ]
    root = [0.0, 5.0, 0.02]
    controls, parents, splines = [root], [-1], []
    for endpoint in ([-4.0, 5.0, 0.02], [4.0, 5.0, 0.02], [0.0, 9.0, 0.02]):
        j = len(controls)
        controls.extend(
            [
                [root[k] + (endpoint[k] - root[k]) / 3 for k in range(3)],
                [root[k] + 2 * (endpoint[k] - root[k]) / 3 for k in range(3)],
                endpoint,
            ]
        )
        parents.extend([0, j + 2, -1])
        splines.append([0, j, j + 1, j + 2])
    identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    return {
        "schema": "curvenet.bundle/1",
        "space": "mesh-local",
        "projection": {"points": points, "faces": faces},
        "network": {
            "schema": 1,
            "points": controls,
            "parents": parents,
            "splines": splines,
        },
        "frames": {"positions": controls, "rotations": [identity for _ in controls]},
    }


def create_demo():
    """Create a bound mesh and controls, and leave an immediately visible bulge.

    This does not clear the scene. Set the selected central knot's Translate Z
    to zero to restore neutral, or move it to change the deformation.
    """
    from maya import cmds
    from maya.api import OpenMaya as om

    data = demo_bundle()
    result = create_rig(data)
    central = result["controls"][0]
    cmds.setAttr(central + ".translateZ", 2.0)
    shape = cmds.listRelatives(
        result["mesh"], shapes=True, noIntermediate=True, fullPath=True
    )[0]
    selection = om.MSelectionList()
    selection.add(shape)
    actual = om.MFnMesh(selection.getDagPath(0)).getPoints(om.MSpace.kObject)
    moved = max(
        (point - om.MPoint(*rest)).length()
        for point, rest in zip(actual, data["projection"]["points"])
    )
    if moved < 1e-5:
        raise RuntimeError(
            "cnProfileMover was created but the surface did not move. Run curvenet_maya.diagnose(%r) and share the output."
            % result["deformer"]
        )
    cmds.select(result["mesh"])
    if not cmds.about(batch=True):
        try:
            cmds.viewFit()
        except RuntimeError:
            pass
    cmds.select(central)
    print(
        "CurveNet demo: node=%s, mesh=%s, active control=%s, max displacement=%g"
        % (result["deformer"], result["mesh"], central, moved)
    )
    result["max_displacement"] = moved
    return result


def diagnose(node):
    """Report node configuration without changing the scene."""
    from maya import cmds
    from maya.api import OpenMaya as om

    if not cmds.objExists(node) or cmds.nodeType(node) != "cnProfileMover":
        raise ValueError("Expected an existing cnProfileMover node")
    selected = om.MSelectionList()
    selected.add(node)
    fn = om.MFnDependencyNode(selected.getDependNode(0))

    def count(attribute, data_type):
        value = fn.findPlug(attribute, False).asMObject()
        return 0 if value.isNull() else len(data_type(value).array())

    report = {
        "node": node,
        "maya_api": cmds.about(apiVersion=True),
        "envelope": cmds.getAttr(node + ".envelope"),
        "geometry": cmds.deformer(node, query=True, geometry=True) or [],
        "projection_configured": not fn.findPlug("projectionMesh", False)
        .asMObject()
        .isNull(),
        "control_points": count("restControlPoints", om.MFnPointArrayData),
        "spline_indices": count("splineIndices", om.MFnIntArrayData),
        "parent_indices": count("controlParents", om.MFnIntArrayData),
        "connected_control_inputs": len(
            cmds.listConnections(
                node + ".controlOffsets", source=True, destination=False, plugs=True
            )
            or []
        ),
    }
    print(json.dumps(report, indent=2))
    return report
