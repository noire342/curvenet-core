"""Run with mayapy tests/maya2023_smoke.py after module registration.

Creates a disposable scene; does not touch a user's current GUI scene.
Writes evidence only after load, native evaluation and save/reopen pass.
"""

import json
import argparse
import sys
import tempfile
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--plugin", help="absolute path to the Maya 2023 .mll")
args = parser.parse_args()
MODULE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_DIR / "scripts"))
import maya.standalone

maya.standalone.initialize(name="python")
from maya import cmds
from maya.api import OpenMaya as om
import curvenet_maya as cn

assert cmds.about(apiVersion=True) // 10000 == 2023, cmds.about(apiVersion=True)
cn.load_plugin(args.plugin)
cmds.file(new=True, force=True)
demo = cn.create_demo()
assert demo["max_displacement"] > 1e-5, "Synthetic demo did not deform"
cmds.file(new=True, force=True)
test_mesh = cmds.polyPlane(
    width=10, height=10, subdivisionsX=10, subdivisionsY=10, axis=(0, 0, 1)
)[0]
test_curves = [
    cmds.curve(degree=1, point=[(-4, 0, 0.02), (0, 0, 0.02)], knot=[0, 1]),
    cmds.curve(degree=1, point=[(0, 0, 0.02), (4, 0, 0.02)], knot=[0, 1]),
    cmds.curve(degree=1, point=[(0, 0, 0.02), (0, 4, 0.02)], knot=[0, 1]),
]
multi = cn.bind_curves(test_mesh, test_curves)
test_shape = cmds.listRelatives(
    test_mesh, shapes=True, noIntermediate=True, fullPath=True
)[0]
selection = om.MSelectionList()
selection.add(test_shape)
test_path = selection.getDagPath(0)
before = om.MFnMesh(test_path).getPoints(om.MSpace.kObject)
cmds.setAttr(multi["controls"][1] + ".translateZ", 1.0)
after = om.MFnMesh(test_path).getPoints(om.MSpace.kObject)
assert max((a - b).length() for a, b in zip(after, before)) > 1e-5, (
    "Multiple curve rig did not deform"
)
curve_report = cn.diagnose(multi["deformer"])
assert curve_report["spline_indices"] == 12, curve_report
cmds.file(new=True, force=True)
n = 6
points = [[x / n, y / n, 0.0] for y in range(n + 1) for x in range(n + 1)]
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
a, b = [0.17, 0.12, 0.02], [0.83, 0.88, 0.02]
cp = [
    a,
    [a[i] + (b[i] - a[i]) / 3 for i in range(3)],
    [a[i] + 2 * (b[i] - a[i]) / 3 for i in range(3)],
    b,
]
identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
bundle = {
    "schema": "curvenet.bundle/1",
    "space": "mesh-local",
    "projection": {"points": points, "faces": faces},
    "network": {
        "schema": 1,
        "points": cp,
        "splines": [[0, 1, 2, 3]],
        "parents": [-1, 0, 3, -1],
    },
    "frames": {"positions": cp, "rotations": [identity] * 4},
}


def positions(mesh):
    shape = cmds.listRelatives(mesh, shapes=True, noIntermediate=True, fullPath=True)[0]
    selection = om.MSelectionList()
    selection.add(shape)
    return om.MFnMesh(selection.getDagPath(0)).getPoints(om.MSpace.kObject)


def error(actual, expected):
    return max((p - om.MPoint(*q)).length() for p, q in zip(actual, expected))


with tempfile.TemporaryDirectory(prefix="curvenet-maya2023-") as directory:
    path = Path(directory) / "grid.json"
    path.write_text(json.dumps(bundle))
    marker = cmds.createNode("transform", name="PreviousUserAction")
    cmds.setAttr(marker + ".tx", 0.4321)
    try:
        cn.import_bundle(path, mesh=marker)
    except ValueError:
        pass
    else:
        raise AssertionError("Accepted a transform without a polygon mesh")
    assert abs(cmds.getAttr(marker + ".tx") - 0.4321) < 1e-10, (
        "Failed import undid a previous user action"
    )
    result = cn.import_bundle(path)
    mesh, node = result["mesh"], result["deformer"]
    assert error(positions(mesh), points) < 1e-7, "Neutral drift"
    for i in (0, 3):
        cmds.setAttr(result["controls"][i] + ".tz", 0.2)
    expected = [[p[0], p[1], p[2] + 0.2] for p in points]
    assert error(positions(mesh), expected) < 1e-6, "Maya row/column channel mapping"
    cmds.setAttr(node + ".envelope", 0.5)
    assert error(positions(mesh), [[p[0], p[1], p[2] + 0.1] for p in points]) < 1e-6, (
        "Envelope"
    )
    cmds.setAttr(node + ".envelope", 0)
    assert error(positions(mesh), points) < 1e-7, "Envelope zero"
    cmds.setAttr(node + ".envelope", 1)
    for i in (0, 3):
        cmds.setAttr(result["controls"][i] + ".tz", 0)
    target = cmds.createNode("transform", name="UpstreamTarget")
    selection = om.MSelectionList()
    selection.add(target)
    om.MFnMesh().create(
        [om.MPoint(p[0], p[1], p[2] + 0.1) for p in points],
        [len(f) for f in faces],
        [i for f in faces for i in f],
        parent=selection.getDependNode(0),
    )
    blend = cmds.blendShape(target, mesh, frontOfChain=True)[0]
    cmds.setAttr(blend + ".weight[0]", 1)
    assert error(positions(mesh), [[p[0], p[1], p[2] + 0.1] for p in points]) < 1e-6, (
        "Upstream neutral transport"
    )
    for i in (0, 3):
        cmds.setAttr(result["controls"][i] + ".tz", 0.2)
    expected = [[p[0], p[1], p[2] + 0.3] for p in points]
    assert error(positions(mesh), expected) < 1e-6, (
        "Correction after upstream blendshape"
    )
    scene = str(Path(directory) / "smoke.ma")
    cmds.file(rename=scene)
    cmds.file(save=True, type="mayaAscii", force=True)
    cmds.file(scene, open=True, force=True)
    assert error(positions(mesh), expected) < 1e-6, "Scene reopen"
    # Descriptor changes must fail cleanly, never index an old cache.
    fn = om.MFnDependencyNode(om.MSelectionList().add(node).getDependNode(0))
    diagnostics = []

    def observe(message, message_type, _client):
        if message_type == om.MCommandMessage.kError:
            diagnostics.append(str(message))

    callback = om.MCommandMessage.addCommandOutputCallback(observe)
    try:
        fn.findPlug("controlParents", False).setMObject(
            om.MFnIntArrayData().create([-1, 999, 3, -1])
        )
        try:
            positions(mesh)
        except RuntimeError:
            pass  # Evaluation failure is permitted only with the expected diagnostic.
    finally:
        om.MMessage.removeCallback(callback)
    assert any(
        "root knot" in message and "CurveNet" in message for message in diagnostics
    ), ("Invalid parent was not explicitly rejected", diagnostics)
    fn.findPlug("controlParents", False).setMObject(
        om.MFnIntArrayData().create([-1, 0, 3, -1])
    )
    assert error(positions(mesh), expected) < 1e-6, "Recovery after invalid descriptor"
report = {
    "maya": cmds.about(version=True),
    "api": cmds.about(apiVersion=True),
    "status": "PASS",
    "checks": [
        "load",
        "visible synthetic demo",
        "three NURBS curve inputs",
        "neutral",
        "direct translation",
        "envelope",
        "upstream blendshape + correction",
        "save/reopen",
        "invalid parent recovery",
    ],
}
report_directory = Path.cwd() / "build"
report_directory.mkdir(exist_ok=True)
(report_directory / "maya2023-runtime-validation.json").write_text(
    json.dumps(report, indent=2)
)
print("MAYA2023_RUNTIME_PASS", report)
maya.standalone.uninitialize()
