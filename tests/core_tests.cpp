#include "curvenet/core.h"
#include "curvenet/rig.h"
#include <Eigen/Geometry>
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>
using namespace cn;
namespace {
Mesh grid(int nx = 6, int ny = 6, bool triangles = false) {
  Mesh m;
  for (int y = 0; y <= ny; ++y)
    for (int x = 0; x <= nx; ++x)
      m.points.emplace_back(double(x) / nx, double(y) / ny, 0);
  for (int y = 0; y < ny; ++y)
    for (int x = 0; x < nx; ++x) {
      int a = y * (nx + 1) + x, b = a + 1, c = b + nx + 1, d = a + nx + 1;
      if (triangles) {
        m.faces.push_back({a, b, c});
        m.faces.push_back({a, c, d});
      } else
        m.faces.push_back({a, b, c, d});
    }
  return m;
}
Net line(Vec3 a = Vec3(.17, .12, .02), Vec3 b = Vec3(.83, .88, .02)) {
  return {{a, a + (b - a) / 3, a + 2 * (b - a) / 3, b}, {{0, 1, 2, 3}}};
}
Net tee() {
  Net n;
  n.points = {Vec3(.5, .5, .01)};
  for (Vec3 p : {Vec3(.13, .5, .01), Vec3(.87, .5, .01), Vec3(.5, .88, .01)}) {
    int j = int(n.points.size());
    Vec3 c = n.points[0];
    n.points.push_back(c + (p - c) / 3);
    n.points.push_back(c + 2 * (p - c) / 3);
    n.points.push_back(p);
    n.splines.push_back({0, j, j + 1, j + 2});
  }
  return n;
}
void check(bool c, const char *msg) {
  if (!c)
    throw std::runtime_error(msg);
}
double error(const std::vector<Vec3> &a, const std::vector<Vec3> &b) {
  double e = 0;
  check(a.size() == b.size(), "point count");
  for (size_t i = 0; i < a.size(); ++i)
    e = std::max(e, (a[i] - b[i]).norm());
  return e;
}
void expectThrow(std::function<void()> f) {
  bool caught = false;
  try {
    f();
  } catch (const std::exception &) {
    caught = true;
  }
  check(caught, "Expected validation failure");
}
} // namespace
int main() {
  int failures = 0;
  auto test = [&](const char *name, std::function<void()> f) {
    try {
      f();
      std::cout << "PASS " << name << '\n';
    } catch (const std::exception &e) {
      ++failures;
      std::cerr << "FAIL " << name << ": " << e.what() << '\n';
    }
  };
  test("polygon Laplacian: constants, symmetry, planar affine precision", [] {
    std::vector<Vec3> p{{0, 0, 0}, {2, 0, 0}, {2, 1, 0}, {0, 1, 0}};
    auto L = polygonLaplacian(p);
    check((L - L.transpose()).norm() < 1e-12, "symmetric");
    check((L * Eigen::VectorXd::Ones(4)).norm() < 1e-12, "constant nullspace");
    Eigen::VectorXd x(4);
    x << 0, 2, 2, 0;
    check(std::abs(x.dot(L * x) - 2) < 1e-12,
          "Dirichlet energy of linear field");
  });
  for (bool tri : {false, true})
    test(tri ? "tri mesh rest identity" : "quad mesh rest identity", [tri] {
      auto m = grid(8, 8, tri);
      auto n = line();
      Solver s;
      s.bind(m, n);
      auto x = s.solve(m.points, n.points);
      check(error(x, m.points) < 1e-8, "Rest pose drift");
      check(s.stats().cutFaces > int(m.faces.size()),
            "Expected actual mesh cuts");
    });
  test("off-surface residual + rigid translation", [] {
    auto m = grid();
    auto n = line();
    Solver s;
    s.bind(m, n);
    Vec3 t(.3, -.15, .7);
    auto q = n.points;
    for (auto &p : q)
      p += t;
    auto expected = m.points;
    for (auto &p : expected)
      p += t;
    check(error(s.solve(m.points, q), expected) < 1e-8,
          "Translation/residual mismatch");
  });
  test("T junction separate side frames + rigid transform", [] {
    auto m = grid(8, 8);
    auto n = tee();
    Solver s;
    s.bind(m, n);
    Mat3 R =
        Eigen::AngleAxisd(.65, Vec3(1, 2, 3).normalized()).toRotationMatrix();
    auto q = n.points;
    auto expected = m.points;
    for (auto &v : q)
      v = R * v + Vec3(.1, .2, .3);
    for (auto &v : expected)
      v = R * v + Vec3(.1, .2, .3);
    check(error(s.solve(m.points, q), expected) < 2e-7,
          "Junction rigid equivariance");
  });
  test("curve coincident with mesh vertices and edges", [] {
    auto m = grid();
    auto n = line(Vec3(0, .5, 0), Vec3(1, .5, 0));
    Solver s;
    s.bind(m, n);
    check(error(s.solve(m.points, n.points), m.points) < 1e-8,
          "Coincident curve drift");
  });
  test("open-ended cut cracks", [] {
    auto m = grid(1, 1);
    auto n = line(Vec3(0, .4, 0), Vec3(.7, .4, 0));
    Solver s;
    s.bind(m, n);
    auto p = n.points;
    for (auto &v : p)
      v.z() += .2;
    auto x = s.solve(m.points, p);
    check(s.stats().corners > 4, "Crack corners were not duplicated");
    for (size_t i = 0; i < x.size(); ++i) {
      if (std::abs(x[i].z() - .2) >= 1e-8)
        std::cout << "  crack vertex " << i << " z=" << x[i].z() << '\n';
      check(std::abs(x[i].z() - .2) < 1e-8, "Crack constraint lost");
    }
  });
  test("isolated sub-face network is removed", [] {
    auto m = grid(1, 1);
    auto n = line(Vec3(.2, .4, 0), Vec3(.7, .4, 0));
    Solver s;
    s.bind(m, n);
    auto p = n.points;
    for (auto &v : p)
      v.z() += 2;
    check(error(s.solve(m.points, p), m.points) < 1e-8,
          "Unresolvable island affected mesh");
    check(s.stats().removedIslands > 0, "Island diagnostics");
  });
  test("upstream curved surface + zero correction", [] {
    auto m = grid(8, 8);
    auto n = tee();
    Solver s;
    s.bind(m, n);
    auto rest = m.points;
    for (auto &v : rest)
      v.z() = .12 * std::sin(v.x() * 4) * std::sin(v.y() * 3);
    auto q = s.warpedControls(rest);
    check(error(s.solve(rest, q), rest) < 1e-7,
          "Upstream pose was modified with zero correction");
  });
  test("repeat solve / deterministic binding cache", [] {
    auto m = grid();
    auto n = tee();
    Solver s;
    s.bind(m, n);
    auto p = n.points;
    p.back().z() += .15;
    auto a = s.solve(m.points, p), b = s.solve(m.points, p);
    check(error(a, b) < 1e-12, "Non-deterministic solve");
    check(s.stats().relativeResidual < 1e-8, "Linear residual");
  });
  test("reject non-manifold, invalid winding, non-finite controls", [] {
    auto m = grid(1, 1);
    auto n = line();
    auto bad = m;
    bad.faces.push_back({1, 0, 3});
    bad.faces.push_back({0, 1, 2});
    expectThrow([&] {
      Solver s;
      s.bind(bad, n);
    });
    bad = m;
    bad.faces.push_back({0, 1, 2});
    expectThrow([&] {
      Solver s;
      s.bind(bad, n);
    });
    n.points[0].x() = NAN;
    expectThrow([&] {
      Solver s;
      s.bind(m, n);
    });
  });
  test("rig binding: descriptor edits invalidate the solver", [] {
    RigSession session;
    auto m = grid(6, 6);
    auto n = line();
    std::vector<int> parents{-1, 0, 3, -1};
    check(session.ensureBinding(m, n, parents, Options{}),
          "Initial bind missing");
    check(!session.ensureBinding(m, n, parents, Options{}),
          "Unchanged descriptor rebound");
    n.points[1].x() += .01;
    check(session.ensureBinding(m, n, parents, Options{}),
          "Control edit reused stale solver");
    m.points[0].z() += .005;
    check(session.ensureBinding(m, n, parents, Options{}),
          "Projection edit reused stale solver");
    check(!session.ensureBinding(m, n, parents, Options{}),
          "Updated descriptor not cached");
    auto invalid = parents;
    invalid[2] = 0;
    expectThrow([&] { session.ensureBinding(m, n, invalid, Options{}); });
    check(!session.ensureBinding(m, n, parents, Options{}),
          "Failed bind destroyed good state");
    n.points.push_back(Vec3(0, 0, 0));
    expectThrow([&] { session.ensureBinding(m, n, parents, Options{}); });
    auto bad = n;
    bad.splines[0][1] = 1000;
    expectThrow([&] { session.ensureBinding(m, bad, parents, Options{}); });
  });
  test("rig evaluation: direct channels, upstream neutral, malformed data", [] {
    RigSession session;
    auto m = grid();
    auto n = line();
    std::vector<ControlOffset> offsets(n.points.size());
    expectThrow([&] { session.evaluate(m, offsets); });
    session.ensureBinding(m, n, {-1, 0, 3, -1}, Options{});
    auto upstream = m;
    for (auto &p : upstream.points)
      p.z() += .03 * p.x();
    check(error(session.evaluate(upstream, offsets), upstream.points) == 0,
          "Neutral rig drift");
    offsets[0].translation.z() = .1;
    offsets[3].translation.z() = .1;
    check(error(session.evaluate(m, offsets), m.points) > .09,
          "Direct controls had no effect");
    offsets[0].translation.x() = NAN;
    expectThrow([&] { session.evaluate(m, offsets); });
    offsets.resize(1);
    expectThrow([&] { session.evaluate(m, offsets); });
    auto bad = m;
    bad.faces[0].pop_back();
    offsets.resize(n.points.size());
    expectThrow([&] { session.evaluate(bad, offsets); });
  });
  test("control rotations reject unbound, wrong size and non-finite rest", [] {
    Solver solver;
    auto m = grid();
    expectThrow([&] { solver.controlRotations(m.points); });
    solver.bind(m, line());
    auto rest = m.points;
    rest.pop_back();
    expectThrow([&] { solver.controlRotations(rest); });
    rest = m.points;
    rest[0].x() = NAN;
    expectThrow([&] { solver.controlRotations(rest); });
  });
  test("rig local frames reproduce known rigid transforms and handle channels",
       [] {
         auto mesh = grid(8, 8);
         auto net = tee();
         std::vector<int> parents{-1, 0, 3, -1, 0, 6, -1, 0, 9, -1};
         RigSession rig;
         rig.ensureBinding(mesh, net, parents, Options{});
         Mat3 rotation = Eigen::AngleAxisd(.6, Vec3(1, 2, 3).normalized())
                             .toRotationMatrix();
         Vec3 translation(.2, -.1, .4);
         auto frame = rig.solver().controlRotations(mesh.points);
         std::vector<ControlOffset> offsets(net.points.size());
         for (int i : {0, 3, 6, 9}) {
           offsets[i].translation =
               frame[i].transpose() *
               (rotation * net.points[i] + translation - net.points[i]);
           offsets[i].linear = frame[i].transpose() * rotation * frame[i];
         }
         auto expected = mesh.points;
         for (auto &p : expected)
           p = rotation * p + translation;
         check(error(rig.evaluate(mesh, offsets), expected) < 1e-7,
               "Local TRS axis or transpose error");
         auto upstream = mesh;
         for (auto &p : upstream.points)
           p = rotation * p + translation;
         offsets.assign(net.points.size(), ControlOffset{});
         auto rotatedFrames = rig.solver().controlRotations(upstream.points);
         Vec3 local(.03, -.04, .05);
         for (int i : {0, 3, 6, 9})
           offsets[i].translation = local;
         expected = upstream.points;
         Vec3 worldDelta = rotatedFrames[0] * local;
         for (auto &p : expected)
           p += worldDelta;
         check(error(rig.evaluate(upstream, offsets), expected) < 1e-7,
               "Upstream local frame mapping error");
         offsets.assign(net.points.size(), ControlOffset{});
         offsets[1].linear = rotation;
         expectThrow([&] { rig.evaluate(mesh, offsets); });
         offsets.assign(net.points.size(), ControlOffset{});
         offsets[0].linear = Mat3::Identity() * 1.2;
         offsets[1].translation = local;
         auto cp = net.points;
         for (int i : {1, 4, 7})
           cp[i] = cp[0] + 1.2 * (cp[i] - cp[0]);
         cp[1] += frame[0] * (1.2 * local);
         Solver reference;
         reference.bind(mesh, net);
         check(error(rig.evaluate(mesh, offsets),
                     reference.solve(mesh.points, cp)) < 1e-9,
               "Parent scale and handle translation composition error");
       });
  std::cout << "Failures: " << failures << '\n';
  return failures ? 1 : 0;
}
