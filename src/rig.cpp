#include "curvenet/rig.h"
#include <stdexcept>

namespace cn {
namespace {
bool samePoints(const std::vector<Vec3> &a, const std::vector<Vec3> &b) {
  if (a.size() != b.size())
    return false;
  for (size_t i = 0; i < a.size(); ++i)
    if (a[i] != b[i])
      return false;
  return true;
}
void validateParents(const Net &net, const std::vector<int> &parents) {
  if (parents.size() != net.points.size())
    throw std::runtime_error("Control-parent count does not match the network");
  for (size_t i = 0; i < parents.size(); ++i) {
    int parent = parents[i];
    if (parent < -1 || parent >= int(parents.size()) ||
        (parent >= 0 && (parent == int(i) || parents[parent] != -1)))
      throw std::runtime_error("A tangent must refer to one root knot");
  }
  for (const auto &s : net.splines)
    if (parents[s[0]] != -1 || parents[s[3]] != -1 || parents[s[1]] != s[0] ||
        parents[s[2]] != s[3])
      throw std::runtime_error(
          "Spline handles must belong to their endpoint knots");
}
} // namespace

bool RigSession::ensureBinding(const Mesh &projection, const Net &network,
                               const std::vector<int> &parents,
                               const Options &options, int revision) {
  if (solver_ && revision == revision_ &&
      options.samplesPerEdge == options_.samplesPerEdge &&
      options.maxSamples == options_.maxSamples && parents == parents_ &&
      network.splines == network_.splines &&
      samePoints(network.points, network_.points) &&
      projection.faces == projection_.faces &&
      samePoints(projection.points, projection_.points))
    return false;
  validate(projection,
           network); // Check spline indices before indexing parents.
  validateParents(network, parents);
  auto candidate = std::make_unique<Solver>();
  candidate->bind(projection, network, options);
  // Finish allocations before committing the descriptor and factorization.
  auto nextProjection = projection;
  auto nextNetwork = network;
  auto nextParents = parents;
  projection_ = std::move(nextProjection);
  network_ = std::move(nextNetwork);
  parents_ = std::move(nextParents);
  options_ = options;
  revision_ = revision;
  solver_ = std::move(candidate);
  return true;
}

std::vector<Vec3>
RigSession::evaluate(const Mesh &upstream,
                     const std::vector<ControlOffset> &offsets) {
  if (!solver_)
    throw std::runtime_error("CurveNet has no valid binding");
  if (upstream.faces != projection_.faces ||
      upstream.points.size() != projection_.points.size())
    throw std::runtime_error("Input topology changed; rebind the network");
  if (offsets.size() != network_.points.size())
    throw std::runtime_error(
        "Animator channels do not match the bound network");
  bool neutral = true;
  for (size_t i = 0; i < offsets.size(); ++i) {
    const auto &offset = offsets[i];
    if (!offset.linear.allFinite() || !offset.translation.allFinite())
      throw std::runtime_error("Animator channels contain non-finite values");
    if (parents_[i] >= 0 && !offset.linear.isApprox(Mat3::Identity(), 1e-12))
      throw std::runtime_error("Tangent handles support translation only");
    neutral = neutral && offset.linear.isApprox(Mat3::Identity(), 1e-12) &&
              offset.translation.norm() < 1e-12;
  }
  for (const auto &point : upstream.points)
    if (!point.allFinite())
      throw std::runtime_error("Upstream surface contains non-finite values");
  if (neutral)
    return upstream.points;
  auto base = solver_->warpedControls(upstream.points);
  auto rotations = solver_->controlRotations(upstream.points);
  auto posed = base;
  for (size_t i = 0; i < parents_.size(); ++i) {
    int p = parents_[i];
    if (p < 0)
      posed[i] = base[i] + rotations[i] * offsets[i].translation;
    else {
      Vec3 arm = rotations[p].transpose() * (base[i] - base[p]);
      posed[i] = base[p] + rotations[p] * (offsets[p].linear *
                                               (arm + offsets[i].translation) +
                                           offsets[p].translation);
    }
  }
  return solver_->solve(upstream.points, posed);
}

const Solver &RigSession::solver() const {
  if (!solver_)
    throw std::runtime_error("CurveNet has no valid binding");
  return *solver_;
}
} // namespace cn
