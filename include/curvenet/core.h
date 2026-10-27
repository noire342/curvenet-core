#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <array>
#include <memory>
#include <string>
#include <vector>

namespace cn {
using Vec3 = Eigen::Vector3d;
using Mat3 = Eigen::Matrix3d;
struct Mesh {
  std::vector<Vec3> points;
  std::vector<std::vector<int>> faces;
};
// Cubic Bezier splines index a shared pool of control points. Endpoints may be
// shared at arbitrary-valence intersections; tangent handles remain per spline.
struct Net {
  std::vector<Vec3> points;
  std::vector<std::array<int, 4>> splines;
};
struct Attachment {
  std::array<int, 3> vertices{};
  Vec3 bary = Vec3::Zero();
  int face = -1;
  Vec3 projected = Vec3::Zero();
};
struct Options {
  double samplesPerEdge = 5.0;
  int maxSamples = 200000;
};
struct Statistics {
  int samples = 0, cutFaces = 0, corners = 0, unknowns = 0,
      pinnedComponents = 0;
  int removedIslands = 0;
  double bindingMs = 0, solveMs = 0, relativeResidual = 0;
};
// Projection and rest poses deliberately remain distinct (paper section 5).
class Solver {
public:
  Solver();
  ~Solver();
  Solver(Solver &&) noexcept;
  Solver &operator=(Solver &&) noexcept;
  void bind(const Mesh &projection, const Net &net,
            const Options &options = {});
  std::vector<Vec3> solve(const std::vector<Vec3> &restPoints,
                          const std::vector<Vec3> &posedControlPoints);
  std::vector<Vec3> warpedControls(const std::vector<Vec3> &restPoints) const;
  std::vector<Mat3> controlRotations(const std::vector<Vec3> &restPoints) const;
  const std::vector<Attachment> &controlAttachments() const;
  const Statistics &stats() const;

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};
Attachment project(const Mesh &, const Vec3 &);
Mat3 attachmentRotation(const Mesh &, const std::vector<Vec3> &,
                        const Attachment &);
Eigen::MatrixXd polygonLaplacian(const std::vector<Vec3> &);
Vec3 bezier(const std::array<Vec3, 4> &, double t);
void validate(const Mesh &, const Net &);
} // namespace cn
