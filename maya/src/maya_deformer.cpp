// Host adapter only; binding and animator-space mapping live in RigSession.
#include "curvenet/rig.h"
#include <cmath>
#include <map>
#include <maya/MArrayDataHandle.h>
#include <maya/MDataBlock.h>
#include <maya/MFnIntArrayData.h>
#include <maya/MFnMatrixAttribute.h>
#include <maya/MFnMesh.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnPlugin.h>
#include <maya/MFnPointArrayData.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MGlobal.h>
#include <maya/MIntArray.h>
#include <maya/MItGeometry.h>
#include <maya/MMatrix.h>
#include <maya/MPointArray.h>
#include <maya/MPxDeformerNode.h>
#include <mutex>
#include <stdexcept>

namespace {
void checked(MStatus s, const char *operation) {
  if (!s)
    throw std::runtime_error(std::string(operation) + ": " +
                             s.errorString().asChar());
}
cn::Mesh readMesh(MObject object) {
  MStatus s;
  MFnMesh fn(object, &s);
  checked(s, "Expected polygon mesh");
  MPointArray points;
  checked(fn.getPoints(points, MSpace::kObject), "Read points");
  MIntArray counts, indices;
  checked(fn.getVertices(counts, indices), "Read topology");
  cn::Mesh out;
  for (unsigned i = 0; i < points.length(); ++i)
    out.points.emplace_back(points[i].x, points[i].y, points[i].z);
  unsigned offset = 0;
  for (unsigned i = 0; i < counts.length(); ++i) {
    std::vector<int> face;
    for (int j = 0; j < counts[i]; ++j)
      face.push_back(indices[offset++]);
    out.faces.push_back(std::move(face));
  }
  return out;
}
MObject inputObject(MDataBlock &data, MObject attribute) {
  MStatus s;
  auto handle = data.inputValue(attribute, &s);
  checked(s, "Read input attribute");
  return handle.data();
}
class CurveNetDeformer final : public MPxDeformerNode {
public:
  static MObject projectionMesh, restControlPoints, splineIndices,
      controlParents, controlOffsets, sampling, bindVersion;
  static MTypeId id;
  static void *creator() { return new CurveNetDeformer; }
  SchedulingType schedulingType() const override { return kSerial; }
  static MStatus initialize() {
    try {
      MFnTypedAttribute t;
      MFnNumericAttribute n;
      MFnMatrixAttribute m;
      MStatus s;
      auto typed = [&](const char *name, const char *shortName,
                       MFnData::Type type) {
        MObject attr = t.create(name, shortName, type, MObject::kNullObj, &s);
        checked(s, name);
        checked(addAttribute(attr), name);
        return attr;
      };
      projectionMesh = typed("projectionMesh", "prj", MFnData::kMesh);
      restControlPoints =
          typed("restControlPoints", "rcp", MFnData::kPointArray);
      splineIndices = typed("splineIndices", "spi", MFnData::kIntArray);
      controlParents = typed("controlParents", "cpa", MFnData::kIntArray);
      controlOffsets =
          m.create("controlOffsets", "cof", MFnMatrixAttribute::kDouble, &s);
      checked(s, "Create control offsets");
      checked(m.setArray(true), "Set array");
      checked(m.setUsesArrayDataBuilder(true), "Set array builder");
      checked(addAttribute(controlOffsets), "Add offsets");
      sampling =
          n.create("samplesPerEdge", "spe", MFnNumericData::kDouble, 5., &s);
      checked(s, "Create sampling");
      checked(n.setMin(1.), "Sampling minimum");
      checked(n.setMax(20.), "Sampling maximum");
      checked(addAttribute(sampling), "Add sampling");
      bindVersion = n.create("bindVersion", "bvr", MFnNumericData::kInt, 0, &s);
      checked(s, "Create bind revision");
      checked(addAttribute(bindVersion), "Add bind revision");
      for (auto attr : {projectionMesh, restControlPoints, splineIndices,
                        controlParents, controlOffsets, sampling, bindVersion})
        checked(attributeAffects(attr, outputGeom), "Register dependency");
      return MS::kSuccess;
    } catch (const std::exception &e) {
      MGlobal::displayError(MString("CurveNet initialization: ") + e.what());
      return MS::kFailure;
    }
  }
  MStatus deform(MDataBlock &data, MItGeometry &it, const MMatrix &,
                 unsigned geometryIndex) override {
    std::lock_guard<std::mutex> lock(mutex_);
    try {
      MStatus s;
      float env = data.inputValue(envelope, &s).asFloat();
      checked(s, "Read envelope");
      if (!std::isfinite(env))
        throw std::runtime_error("Envelope is not finite");
      if (env == 0)
        return MS::kSuccess;
      MObject projection = inputObject(data, projectionMesh);
      if (projection.isNull())
        return MS::kSuccess; // Newly created, unbound node.
      MFnPointArrayData cpData(inputObject(data, restControlPoints), &s);
      checked(s, "Read restControlPoints");
      MPointArray cp = cpData.array(&s);
      checked(s, "Read point array");
      MFnIntArrayData ixData(inputObject(data, splineIndices), &s);
      checked(s, "Read splineIndices");
      MIntArray ix = ixData.array(&s);
      checked(s, "Read spline index array");
      MFnIntArrayData parentData(inputObject(data, controlParents), &s);
      checked(s, "Read controlParents");
      MIntArray parent = parentData.array(&s);
      checked(s, "Read parent array");
      if (!cp.length() || !ix.length() || ix.length() % 4)
        throw std::runtime_error("Empty network or incomplete spline indices");
      cn::Net net;
      std::vector<int> parents;
      for (unsigned i = 0; i < cp.length(); ++i)
        net.points.emplace_back(cp[i].x, cp[i].y, cp[i].z);
      for (unsigned i = 0; i < parent.length(); ++i)
        parents.push_back(parent[i]);
      for (unsigned i = 0; i < ix.length(); i += 4)
        net.splines.push_back({ix[i], ix[i + 1], ix[i + 2], ix[i + 3]});
      double density = data.inputValue(sampling, &s).asDouble();
      checked(s, "Read density");
      int revision = data.inputValue(bindVersion, &s).asInt();
      checked(s, "Read binding revision");
      auto &session = sessions_[geometryIndex];
      session.ensureBinding(readMesh(projection), net, parents,
                            cn::Options{density}, revision);
      auto inputs = data.inputArrayValue(input, &s);
      checked(s, "Read geometry inputs");
      checked(inputs.jumpToElement(geometryIndex), "Find input geometry index");
      auto incoming = inputs.inputValue(&s);
      checked(s, "Read input geometry");
      cn::Mesh upstream = readMesh(incoming.child(inputGeom).asMesh());
      std::vector<cn::ControlOffset> offsets(cp.length());
      auto array = data.inputArrayValue(controlOffsets, &s);
      checked(s, "Read animator channels");
      unsigned channelCount = array.elementCount(&s);
      checked(s, "Count animator channels");
      if (channelCount)
        checked(array.jumpToArrayElement(0), "Find first animator channel");
      for (unsigned i = 0; i < channelCount; ++i) {
        unsigned index = array.elementIndex(&s);
        checked(s, "Read animator index");
        if (index >= offsets.size())
          throw std::runtime_error("Animator index exceeds control pool");
        MMatrix matrix = array.inputValue(&s).asMatrix();
        checked(s, "Read animator matrix");
        if (matrix(0, 3) != 0 || matrix(1, 3) != 0 || matrix(2, 3) != 0 ||
            matrix(3, 3) != 1)
          throw std::runtime_error("Animator matrix must be affine");
        for (int row = 0; row < 3; ++row)
          for (int col = 0; col < 3; ++col)
            offsets[index].linear(row, col) = matrix(col, row);
        offsets[index].translation =
            cn::Vec3(matrix(3, 0), matrix(3, 1), matrix(3, 2));
        if (i + 1 < channelCount)
          checked(array.next(), "Advance animator array");
      }
      auto result = session.evaluate(upstream, offsets);
      for (; !it.isDone(); checked(it.next(), "Advance geometry iterator")) {
        int index = it.index(&s);
        checked(s, "Read vertex index");
        if (index < 0 || index >= int(result.size()))
          throw std::runtime_error("Geometry index outside mesh");
        double weight = env * weightValue(data, geometryIndex, index);
        if (!std::isfinite(weight))
          throw std::runtime_error("Deformer weight is not finite");
        if (weight == 0)
          continue;
        cn::Vec3 p = upstream.points[index] +
                     weight * (result[index] - upstream.points[index]);
        checked(it.setPosition(MPoint(p.x(), p.y(), p.z())),
                "Write deformed point");
      }
      lastError_.clear();
      return MS::kSuccess;
    } catch (const std::exception &e) {
      if (lastError_ != e.what()) {
        lastError_ = e.what();
        MGlobal::displayError(MString("CurveNet: ") + e.what());
      }
      return MS::kFailure;
    }
  }

private:
  std::mutex mutex_;
  std::map<unsigned, cn::RigSession> sessions_;
  std::string lastError_;
};
MTypeId CurveNetDeformer::id(0x0007F301);
MObject CurveNetDeformer::projectionMesh, CurveNetDeformer::restControlPoints,
    CurveNetDeformer::splineIndices, CurveNetDeformer::controlParents,
    CurveNetDeformer::controlOffsets, CurveNetDeformer::sampling,
    CurveNetDeformer::bindVersion;
} // namespace
MStatus initializePlugin(MObject object) {
  MFnPlugin p(object, "CurveNet Research", "0.2.1", "Any");
  return p.registerNode("cnProfileMover", CurveNetDeformer::id,
                        CurveNetDeformer::creator, CurveNetDeformer::initialize,
                        MPxNode::kDeformerNode);
}
MStatus uninitializePlugin(MObject object) {
  return MFnPlugin(object).deregisterNode(CurveNetDeformer::id);
}
