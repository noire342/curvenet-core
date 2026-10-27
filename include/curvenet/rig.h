#pragma once
#include "curvenet/core.h"

namespace cn {
// Host-independent animator channels, in each surface-attached knot frame.
struct ControlOffset {
  Mat3 linear = Mat3::Identity();
  Vec3 translation = Vec3::Zero();
};

// Owns one immutable binding. Hosts supply complete descriptors on evaluation;
// changing any binding data cannot accidentally reuse a stale factorization.
class RigSession {
public:
  bool ensureBinding(const Mesh &projection, const Net &network,
                     const std::vector<int> &parents, const Options &options,
                     int revision = 0);
  std::vector<Vec3> evaluate(const Mesh &upstream,
                             const std::vector<ControlOffset> &offsets);
  const Solver &solver() const;

private:
  Mesh projection_;
  Net network_;
  std::vector<int> parents_;
  Options options_;
  int revision_ = 0;
  std::unique_ptr<Solver> solver_;
};
} // namespace cn
