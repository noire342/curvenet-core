#include "curvenet/core.h"
#include <iostream>

int main() {
  cn::Mesh mesh;
  constexpr int n = 6;
  for (int y = 0; y <= n; ++y)
    for (int x = 0; x <= n; ++x)
      mesh.points.emplace_back(double(x) / n, double(y) / n, 0);
  for (int y = 0; y < n; ++y)
    for (int x = 0; x < n; ++x) {
      int a = y * (n + 1) + x;
      mesh.faces.push_back({a, a + 1, a + n + 2, a + n + 1});
    }
  cn::Vec3 a(.17, .12, .02), b(.83, .88, .02);
  cn::Net network{{a, a + (b-a)/3, a + 2*(b-a)/3, b}, {{0,1,2,3}}};
  cn::Solver solver;
  solver.bind(mesh, network);
  auto posed = network.points;
  for (auto &control : posed) control.z() += .2;
  auto result = solver.solve(mesh.points, posed);
  double max_error = 0;
  for (size_t i = 0; i < result.size(); ++i)
    max_error = std::max(max_error, (result[i] - (mesh.points[i] + cn::Vec3(0, 0, .2))).norm());
  std::cout << "vertices: " << result.size() << ", translation error: " << max_error << '\n';
  return max_error < 1e-7 ? 0 : 1;
}
