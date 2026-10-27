#include "curvenet/core.h"
#include <Eigen/Geometry>
#include <Eigen/SVD>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <stdexcept>
#include <unordered_map>

namespace cn {
namespace {
using Clock = std::chrono::steady_clock;
using Sp = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;
constexpr double tiny = 1e-12;
double ms(Clock::time_point a) {
  return std::chrono::duration<double, std::milli>(Clock::now() - a).count();
}
void require(bool ok, const std::string &msg) {
  if (!ok)
    throw std::runtime_error(msg);
}
Vec3 unit(const Vec3 &v) {
  require(v.norm() > tiny,
          "Collapsed curve segment or degenerate surface frame");
  return v.normalized();
}
Mat3 skew(const Vec3 &n) {
  Mat3 m;
  m << 0, -n.z(), n.y(), n.z(), 0, -n.x(), -n.y(), n.x(), 0;
  return m;
}
std::pair<int, int> key(int a, int b) { return std::minmax(a, b); }
Vec3 normal(const Mesh &m, int f) {
  Vec3 n = Vec3::Zero();
  const auto &p = m.faces[f];
  for (size_t i = 0; i < p.size(); ++i)
    n += m.points[p[i]].cross(m.points[p[(i + 1) % p.size()]]);
  return unit(n);
}
double cross2(const Eigen::Vector2d &a, const Eigen::Vector2d &b) {
  return a.x() * b.y() - a.y() * b.x();
}
// Closest triangle point, including edge/vertex regions. The barycentric
// attachment is reused for projection/rest separation; UV seams do not matter.
Vec3 closestTri(const Vec3 &p, const Vec3 &a, const Vec3 &b, const Vec3 &c) {
  Vec3 ab = b - a, ac = c - a, ap = p - a;
  double d1 = ab.dot(ap), d2 = ac.dot(ap);
  if (d1 <= 0 && d2 <= 0)
    return Vec3(1, 0, 0);
  Vec3 bp = p - b;
  double d3 = ab.dot(bp), d4 = ac.dot(bp);
  if (d3 >= 0 && d4 <= d3)
    return Vec3(0, 1, 0);
  double vc = d1 * d4 - d3 * d2;
  if (vc <= 0 && d1 >= 0 && d3 <= 0) {
    double v = d1 / (d1 - d3);
    return Vec3(1 - v, v, 0);
  }
  Vec3 cp = p - c;
  double d5 = ab.dot(cp), d6 = ac.dot(cp);
  if (d6 >= 0 && d5 <= d6)
    return Vec3(0, 0, 1);
  double vb = d5 * d2 - d1 * d6;
  if (vb <= 0 && d2 >= 0 && d6 <= 0) {
    double w = d2 / (d2 - d6);
    return Vec3(1 - w, 0, w);
  }
  double va = d3 * d6 - d5 * d4;
  if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
    double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    return Vec3(0, 1 - w, w);
  }
  double den = 1.0 / (va + vb + vc), v = vb * den, w = vc * den;
  return Vec3(1 - v - w, v, w);
}
struct Tri {
  std::array<int, 3> v;
  int face;
  Vec3 lo, hi, center;
};
struct BVH {
  struct Branch {
    Vec3 lo, hi;
    int begin, end, left = -1, right = -1;
  };
  const Mesh *mesh = nullptr;
  std::vector<Tri> tris;
  std::vector<Branch> tree;
  explicit BVH(const Mesh &m) : mesh(&m) {
    // Ear clipping in the polygon's dominant tangent plane. Numerical
    // integration still uses the original polygon, never its triangles.
    for (int f = 0; f < int(m.faces.size()); ++f) {
      auto ids = m.faces[f];
      Vec3 n = normal(m, f), u = n.unitOrthogonal(), v = n.cross(u);
      auto xy = [&](int i) {
        return Eigen::Vector2d(m.points[i].dot(u), m.points[i].dot(v));
      };
      auto add = [&](int a, int b, int c) {
        Tri t{{a, b, c},
              f,
              m.points[a].cwiseMin(m.points[b]).cwiseMin(m.points[c]),
              m.points[a].cwiseMax(m.points[b]).cwiseMax(m.points[c]),
              (m.points[a] + m.points[b] + m.points[c]) / 3};
        tris.push_back(t);
      };
      int guard = 0;
      while (ids.size() > 3) {
        bool found = false;
        for (int j = 0; j < int(ids.size()); ++j) {
          int a = ids[(j + ids.size() - 1) % ids.size()], b = ids[j],
              c = ids[(j + 1) % ids.size()];
          auto A = xy(a), B = xy(b), C = xy(c);
          if (cross2(B - A, C - B) <= tiny)
            continue;
          bool inside = false;
          for (int k : ids)
            if (k != a && k != b && k != c) {
              auto P = xy(k);
              if (cross2(B - A, P - A) >= -tiny &&
                  cross2(C - B, P - B) >= -tiny &&
                  cross2(A - C, P - C) >= -tiny) {
                inside = true;
                break;
              }
            }
          if (!inside) {
            add(a, b, c);
            ids.erase(ids.begin() + j);
            found = true;
            break;
          }
        }
        require(found && ++guard < 10000,
                "Cannot triangulate polygon " + std::to_string(f) +
                    "; check self-intersections or duplicate vertices");
      }
      add(ids[0], ids[1], ids[2]);
    }
    build(0, int(tris.size()));
  }
  int build(int a, int b) {
    Branch n;
    n.begin = a;
    n.end = b;
    n.lo = Vec3::Constant(std::numeric_limits<double>::infinity());
    n.hi = -n.lo;
    for (int i = a; i < b; ++i) {
      n.lo = n.lo.cwiseMin(tris[i].lo);
      n.hi = n.hi.cwiseMax(tris[i].hi);
    }
    int id = int(tree.size());
    tree.push_back(n);
    if (b - a > 8) {
      Eigen::Index axis;
      (n.hi - n.lo).maxCoeff(&axis);
      int mid = (a + b) / 2;
      std::nth_element(tris.begin() + a, tris.begin() + mid, tris.begin() + b,
                       [&](const Tri &x, const Tri &y) {
                         return x.center[axis] < y.center[axis];
                       });
      int l = build(a, mid), r = build(mid, b);
      tree[id].left = l;
      tree[id].right = r;
    }
    return id;
  }
  Attachment query(const Vec3 &p) const {
    Attachment hit;
    double best = std::numeric_limits<double>::infinity();
    std::vector<int> stack{0};
    while (!stack.empty()) {
      int id = stack.back();
      stack.pop_back();
      const auto &n = tree[id];
      Vec3 d =
          (n.lo - p).cwiseMax(Vec3::Zero()) + (p - n.hi).cwiseMax(Vec3::Zero());
      if (d.squaredNorm() > best)
        continue;
      if (n.left >= 0) {
        stack.push_back(n.left);
        stack.push_back(n.right);
        continue;
      }
      for (int i = n.begin; i < n.end; ++i) {
        const auto &t = tris[i];
        Vec3 bc = closestTri(p, mesh->points[t.v[0]], mesh->points[t.v[1]],
                             mesh->points[t.v[2]]);
        Vec3 q = bc[0] * mesh->points[t.v[0]] + bc[1] * mesh->points[t.v[1]] +
                 bc[2] * mesh->points[t.v[2]];
        double d2 = (p - q).squaredNorm();
        if (d2 < best) {
          best = d2;
          hit = {t.v, bc, t.face, q};
        }
      }
    }
    return hit;
  }
};
Vec3 attached(const std::vector<Vec3> &p, const Attachment &a) {
  return a.bary[0] * p[a.vertices[0]] + a.bary[1] * p[a.vertices[1]] +
         a.bary[2] * p[a.vertices[2]];
}
Mat3 triangleBasis(const std::vector<Vec3> &p, const Attachment &a) {
  Mat3 b;
  Vec3 e = p[a.vertices[1]] - p[a.vertices[0]],
       f = p[a.vertices[2]] - p[a.vertices[0]];
  b.col(0) = e;
  b.col(1) = f;
  b.col(2) = unit(e.cross(f)) * std::sqrt(e.cross(f).norm());
  return b;
}
std::vector<double> arcParameters(const std::array<Vec3, 4> &p, int n) {
  int steps = std::max(64, n * 4);
  std::vector<double> lengths(steps + 1, 0);
  Vec3 prev = p[0];
  for (int i = 1; i <= steps; ++i) {
    Vec3 v = bezier(p, double(i) / steps);
    lengths[i] = lengths[i - 1] + (v - prev).norm();
    prev = v;
  }
  require(lengths.back() > tiny, "Zero length Bezier spline");
  std::vector<double> out(n + 1);
  for (int i = 1; i < n; ++i) {
    double l = lengths.back() * i / n;
    int j = int(std::lower_bound(lengths.begin(), lengths.end(), l) -
                lengths.begin());
    out[i] =
        (j - 1 +
         (l - lengths[j - 1]) / std::max(tiny, lengths[j] - lengths[j - 1])) /
        steps;
  }
  out.back() = 1;
  return out;
}
struct Sample {
  Vec3 q;
  Attachment att;
};
struct Segment {
  int a, b, spline;
};
struct Curve {
  std::vector<int> seg, sign;
  int start, end;
  bool closed = false;
};
struct Piece {
  Vec3 a, b;
  Attachment aa, bb;
  int seg;
  double t0, t1;
};
struct Label {
  int seg = -1;
  double alpha = 0;
  int side = 0;
};
struct CutNode {
  Vec3 p;
  Attachment att;
  int vertex = -1;
  Label fallback;
};
struct Edge {
  int a, b, seg = -1;
  double t0 = 0, t1 = 1;
  bool boundary = false;
};
struct Corner {
  int vertex = -1, unknown = -1;
  Attachment att;
  std::vector<std::pair<int, double>> c;
};
struct Face {
  std::vector<Corner> corners;
  Eigen::MatrixXd L;
};
struct SideFrame {
  Mat3 plus = Mat3::Identity(), minus = Mat3::Identity();
};
} // namespace

Vec3 bezier(const std::array<Vec3, 4> &p, double t) {
  double s = 1 - t;
  return s * s * s * p[0] + 3 * s * s * t * p[1] + 3 * s * t * t * p[2] +
         t * t * t * p[3];
}
Attachment project(const Mesh &m, const Vec3 &p) { return BVH(m).query(p); }
Mat3 attachmentRotation(const Mesh &m, const std::vector<Vec3> &rest,
                        const Attachment &a) {
  Mat3 F = triangleBasis(rest, a) * triangleBasis(m.points, a).inverse();
  Eigen::JacobiSVD<Mat3> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Mat3 d = Mat3::Identity();
  d(2, 2) =
      (svd.matrixU() * svd.matrixV().transpose()).determinant() < 0 ? -1 : 1;
  return svd.matrixU() * d * svd.matrixV().transpose();
}
Eigen::MatrixXd polygonLaplacian(const std::vector<Vec3> &p) {
  int n = int(p.size());
  require(n >= 3, "Polygon needs at least three corners");
  Eigen::MatrixXd X(n, 3), D = Eigen::MatrixXd::Zero(n, n), A = D;
  Vec3 area = Vec3::Zero();
  for (int i = 0; i < n; ++i) {
    X.row(i) = p[i];
    D(i, i) = -1;
    D(i, (i + 1) % n) = 1;
    A(i, i) = A(i, (i + 1) % n) = .5;
    area += .5 * p[i].cross(p[(i + 1) % n]);
  }
  double a = area.norm();
  require(a > tiny, "Zero area cut polygon");
  Eigen::MatrixXd E = D * X, G = (-1 / a) * skew(area / a) * E.transpose() * A,
                  Q = D - E * G;
  return a * G.transpose() * G + Q.transpose() * Q;
}
void validate(const Mesh &m, const Net &net) {
  require(!m.points.empty() && !m.faces.empty(), "Empty projection mesh");
  for (const auto &p : m.points)
    require(p.allFinite(), "Non-finite mesh position");
  std::map<std::pair<int, int>, std::vector<int>> edges;
  for (int f = 0; f < int(m.faces.size()); ++f) {
    const auto &p = m.faces[f];
    require(p.size() >= 3, "Face has fewer than three vertices");
    std::set<int> used;
    for (size_t i = 0; i < p.size(); ++i) {
      int a = p[i], b = p[(i + 1) % p.size()];
      require(a >= 0 && a < int(m.points.size()) && b >= 0 &&
                  b < int(m.points.size()),
              "Invalid mesh vertex index");
      require(used.insert(a).second, "Face contains repeated vertices");
      require((m.points[a] - m.points[b]).norm() > tiny,
              "Zero length mesh edge");
      auto &e = edges[key(a, b)];
      e.push_back(a < b ? 1 : -1);
      require(e.size() <= 2,
              "Non-manifold mesh edge; repair the input mesh first");
    }
    normal(m, f);
  }
  for (const auto &e : edges)
    if (e.second.size() == 2)
      require(e.second[0] != e.second[1], "Inconsistent polygon winding");
  require(!net.splines.empty(), "The curve network is empty");
  for (const auto &p : net.points)
    require(p.allFinite(), "Non-finite curve point");
  std::set<std::pair<int, int>> spans;
  for (const auto &s : net.splines) {
    for (int i : s)
      require(i >= 0 && i < int(net.points.size()),
              "Invalid curve control index");
    require(s[0] != s[3], "A closed curve needs at least two Bezier spans");
    require(spans.insert(key(s[0], s[3])).second,
            "Duplicate spline endpoints; split closed curves into at least "
            "three spans");
  }
}

struct Solver::Impl {
  Mesh mesh;
  Net net;
  Options options;
  Statistics stats;
  std::vector<Sample> samples;
  std::vector<Segment> segments;
  std::vector<std::vector<int>> splineSamples;
  std::vector<Curve> curves;
  std::vector<Attachment> cpAttachments;
  std::map<int, std::vector<int>> junctions;
  std::vector<Face> faces;
  std::vector<int> unknown, freeVertices, pins;
  Sp A, B;
  Eigen::SimplicialLDLT<Sp> factor;
  double eps = 1e-8;

  std::vector<Vec3> warp(const std::vector<Vec3> &rest) const {
    require(rest.size() == mesh.points.size(),
            "Mesh vertex count changed; rebind the network");
    std::vector<Vec3> out(net.points.size());
    for (int i = 0; i < int(out.size()); ++i) {
      const auto &a = cpAttachments[i];
      Mat3 F = triangleBasis(rest, a) * triangleBasis(mesh.points, a).inverse();
      out[i] = attached(rest, a) + F * (net.points[i] - a.projected);
    }
    return out;
  }
  void makeCurves() {
    std::vector<std::vector<int>> adj(net.points.size());
    for (int s = 0; s < int(net.splines.size()); ++s) {
      adj[net.splines[s][0]].push_back(s);
      adj[net.splines[s][3]].push_back(s);
    }
    std::vector<bool> seen(net.splines.size(), false);
    std::vector<int> first(net.splines.size());
    int off = 0;
    for (int i = 0; i < int(first.size()); ++i) {
      first[i] = off;
      off += int(splineSamples[i].size()) - 1;
    }
    auto follow = [&](int start, int initial) {
      Curve c;
      c.start = start;
      int knot = start, s = initial;
      while (!seen[s]) {
        seen[s] = true;
        const auto &sp = net.splines[s];
        int sign = sp[0] == knot ? 1 : -1;
        int count = int(splineSamples[s].size()) - 1;
        for (int j = 0; j < count; ++j) {
          c.seg.push_back(first[s] + (sign == 1 ? j : count - 1 - j));
          c.sign.push_back(sign);
        }
        knot = sign == 1 ? sp[3] : sp[0];
        if (knot == start) {
          c.closed = true;
          break;
        }
        if (adj[knot].size() != 2)
          break;
        s = adj[knot][0] == s ? adj[knot][1] : adj[knot][0];
      }
      c.end = knot;
      if (adj[c.start].size() < 3 && adj[c.end].size() >= 3) {
        std::reverse(c.seg.begin(), c.seg.end());
        std::reverse(c.sign.begin(), c.sign.end());
        for (auto &x : c.sign)
          x = -x;
        std::swap(c.start, c.end);
      }
      curves.push_back(c);
    };
    for (int k = 0; k < int(adj.size()); ++k)
      if (adj[k].size() != 2)
        for (int s : adj[k])
          if (!seen[s])
            follow(k, s);
    for (int s = 0; s < int(seen.size()); ++s)
      if (!seen[s])
        follow(net.splines[s][0], s);
    for (int k = 0; k < int(adj.size()); ++k)
      if (adj[k].size() >= 3) {
        auto list = adj[k];
        Vec3 n = normal(mesh, cpAttachments[k].face), u = n.unitOrthogonal(),
             v = n.cross(u);
        auto direction = [&](int s) {
          const auto &sp = net.splines[s];
          int g = sp[0] == k ? first[s]
                             : first[s] + int(splineSamples[s].size()) - 2;
          const auto &seg = segments[g];
          Vec3 d = samples[seg.b].q - samples[seg.a].q;
          return sp[0] == k ? d : -d;
        };
        std::sort(list.begin(), list.end(), [&](int a, int b) {
          Vec3 da = direction(a), db = direction(b);
          return std::atan2(da.dot(v), da.dot(u)) <
                 std::atan2(db.dot(v), db.dot(u));
        });
        for (int &s : list)
          s = net.splines[s][0] == k
                  ? first[s]
                  : -(first[s] + int(splineSamples[s].size()) - 2) - 1;
        junctions[k] = list;
      }
  }
  std::vector<SideFrame> frames(const std::vector<Vec3> &q) const {
    std::vector<Vec3> tangent(segments.size());
    std::vector<double> len(segments.size());
    for (int i = 0; i < int(segments.size()); ++i) {
      Vec3 d = q[segments[i].b] - q[segments[i].a];
      len[i] = d.norm();
      tangent[i] = unit(d);
    }
    struct Ends {
      Vec3 np, nm;
      double wp, wm;
    };
    std::map<std::pair<int, int>, Ends> end;
    for (const auto &j : junctions) {
      const auto &ids = j.second;
      int n = int(ids.size());
      std::vector<Vec3> t(n), cross(n), norm(n);
      std::vector<double> l(n);
      for (int i = 0; i < n; ++i) {
        int s = ids[i] >= 0 ? ids[i] : -ids[i] - 1;
        t[i] = ids[i] >= 0 ? tangent[s] : -tangent[s];
        l[i] = len[s];
      }
      for (int i = 0; i < n; ++i)
        cross[i] = t[i].cross(t[(i + 1) % n]);
      for (int i = 0; i < n; ++i) {
        Vec3 c = cross[i];
        if (c.norm() < 1e-8)
          c = cross[(i + 1) % n] + cross[(i + n - 1) % n];
        norm[i] = unit(c);
      }
      for (int i = 0; i < n; ++i) {
        int prev = (i + n - 1) % n;
        end[{j.first, ids[i]}] = {norm[i], norm[prev],
                                  l[i] +
                                      cross[i].norm() * (l[(i + 1) % n] - l[i]),
                                  l[i] + cross[prev].norm() * (l[prev] - l[i])};
      }
    }
    std::vector<SideFrame> out(segments.size());
    for (const auto &c : curves) {
      int n = int(c.seg.size());
      bool hasStart = junctions.count(c.start) > 0,
           hasEnd = junctions.count(c.end) > 0;
      if (!hasStart && !hasEnd) {
        for (int s : c.seg) {
          out[s].plus = out[s].minus = Mat3::Zero();
          out[s].plus.col(0) = tangent[s] * len[s];
          out[s].plus.col(1) = tangent[s].unitOrthogonal() * len[s];
          out[s].plus.col(2) = tangent[s].cross(out[s].plus.col(1));
          out[s].minus = out[s].plus;
        }
        continue;
      }
      int startId = c.sign[0] > 0 ? c.seg[0] : -c.seg[0] - 1;
      Ends a = end.at({c.start, startId});
      Ends b = a;
      if (hasEnd) {
        int last = c.seg.back(), dir = -c.sign.back();
        int id = dir > 0 ? last : -last - 1;
        Ends temp = end.at({c.end, id});
        b = {temp.nm, temp.np, temp.wm, temp.wp};
      }
      std::vector<Vec3> ts(n);
      std::vector<double> arclen(n);
      double total = 0;
      for (int i = 0; i < n; ++i) {
        ts[i] = tangent[c.seg[i]] * c.sign[i];
        arclen[i] = total;
        total += len[c.seg[i]];
      }
      for (int side = 0; side < 2; ++side) {
        std::vector<Vec3> normals(n);
        normals[0] = side == 0 ? a.np : a.nm;
        for (int i = 1; i < n; ++i)
          normals[i] = Eigen::Quaterniond::FromTwoVectors(ts[i - 1], ts[i]) *
                       normals[i - 1];
        Vec3 target = side == 0 ? b.np : b.nm;
        double theta =
            hasEnd ? std::atan2(ts.back().dot(normals.back().cross(target)),
                                normals.back().dot(target))
                   : 0;
        for (int i = 0; i < n; ++i) {
          double alpha =
              (n > 1 && arclen.back() > tiny) ? arclen[i] / arclen.back() : 0;
          Vec3 nn = Eigen::AngleAxisd(alpha * theta, ts[i]) * normals[i];
          nn = unit(nn - ts[i] * ts[i].dot(nn));
          double w = (1 - alpha) * (side == 0 ? a.wp : a.wm) +
                     alpha * (side == 0 ? b.wp : b.wm);
          require(w > tiny, "Collapsed curvenet width");
          Mat3 frame;
          frame.col(0) = ts[i] * len[c.seg[i]];
          frame.col(1) = nn.cross(ts[i]) * w;
          frame.col(2) = nn * std::sqrt(len[c.seg[i]] * w);
          int which = c.sign[i] > 0 ? side : 1 - side;
          if (which == 0)
            out[c.seg[i]].plus = frame;
          else
            out[c.seg[i]].minus = frame;
        }
      }
    }
    return out;
  }
  std::vector<SideFrame> gradients(const std::vector<Vec3> &rest,
                                   const std::vector<Vec3> &posed) const {
    auto a = frames(rest), b = frames(posed);
    std::vector<bool> isolated(segments.size(), false);
    for (const auto &c : curves)
      if (!junctions.count(c.start) && !junctions.count(c.end))
        for (int s : c.seg)
          isolated[s] = true;
    for (int s = 0; s < int(a.size()); ++s) {
      if (isolated[s]) {
        Vec3 x = rest[segments[s].b] - rest[segments[s].a],
             y = posed[segments[s].b] - posed[segments[s].a];
        Mat3 F = Eigen::Quaterniond::FromTwoVectors(x, y).toRotationMatrix() *
                 (y.norm() / x.norm());
        b[s].plus = b[s].minus = F;
      } else {
        b[s].plus = b[s].plus * a[s].plus.inverse();
        b[s].minus = b[s].minus * a[s].minus.inverse();
      }
    }
    return b;
  }
  void cutting(const BVH &bvh) {
    std::map<std::pair<int, int>, std::vector<int>> edgeFaces;
    std::vector<std::vector<std::pair<int, std::pair<int, int>>>> neighbors(
        mesh.faces.size());
    std::vector<std::vector<int>> vertexFaces(mesh.points.size());
    for (int f = 0; f < int(mesh.faces.size()); ++f) {
      const auto &p = mesh.faces[f];
      for (size_t j = 0; j < p.size(); ++j) {
        edgeFaces[key(p[j], p[(j + 1) % p.size()])].push_back(f);
        vertexFaces[p[j]].push_back(f);
      }
    }
    for (const auto &e : edgeFaces)
      if (e.second.size() == 2) {
        int a = e.second[0], b = e.second[1];
        neighbors[a].push_back({b, e.first});
        neighbors[b].push_back({a, e.first});
      }
    std::vector<std::vector<Piece>> pieces(mesh.faces.size());
    auto attachment = [&](const Vec3 &p, int face) {
      Attachment a = bvh.query(p);
      if (a.face != face) {
        const auto &ids = mesh.faces[face];
        for (size_t j = 1; j + 1 < ids.size(); ++j) {
          auto bc = closestTri(p, mesh.points[ids[0]], mesh.points[ids[j]],
                               mesh.points[ids[j + 1]]);
          Vec3 v = bc[0] * mesh.points[ids[0]] + bc[1] * mesh.points[ids[j]] +
                   bc[2] * mesh.points[ids[j + 1]];
          if ((p - v).norm() < eps * 4)
            return Attachment{{ids[0], ids[j], ids[j + 1]}, bc, face, v};
        }
      }
      return a;
    };
    auto incident = [&](const Attachment &a) {
      std::vector<int> fs{a.face};
      const auto &p = mesh.faces[a.face];
      for (int id : p)
        if ((a.projected - mesh.points[id]).norm() < eps * .1)
          return vertexFaces[id];
      for (size_t i = 0; i < p.size(); ++i) {
        Vec3 x = mesh.points[p[i]], d = mesh.points[p[(i + 1) % p.size()]] - x;
        double t = (a.projected - x).dot(d) / d.squaredNorm();
        if (t >= -1e-7 && t <= 1.0000001 &&
            (a.projected - (x + std::clamp(t, 0., 1.) * d)).norm() < eps * .1) {
          for (int f : edgeFaces.at(key(p[i], p[(i + 1) % p.size()])))
            if (std::find(fs.begin(), fs.end(), f) == fs.end())
              fs.push_back(f);
        }
      }
      return fs;
    };
    for (int s = 0; s < int(segments.size()); ++s) {
      const auto &seg = segments[s];
      const auto &sa = samples[seg.a];
      const auto &sb = samples[seg.b];
      Vec3 p = sa.att.projected, q = sb.att.projected;
      if ((q - p).norm() < eps * .01)
        continue;
      auto af = incident(sa.att), bf = incident(sb.att);
      int common = -1;
      for (int a : af)
        if (std::find(bf.begin(), bf.end(), a) != bf.end()) {
          common = a;
          break;
        }
      if (common >= 0) {
        pieces[common].push_back({p, q, sa.att, sb.att, s, 0, 1});
        // A curve lying exactly on an original edge constrains both faces.
        for (int f : af)
          if (f != common && std::find(bf.begin(), bf.end(), f) != bf.end())
            pieces[f].push_back({p, q, sa.att, sb.att, s, 0, 1});
        continue;
      }
      // Dense sampling makes paths short. Unfold candidate dual-graph
      // paths and intersect the straight segment with every shared edge.
      std::vector<int> distance(mesh.faces.size(), -1);
      std::queue<int> queue;
      for (int f : bf) {
        distance[f] = 0;
        queue.push(f);
      }
      int minimum = 100000;
      while (!queue.empty()) {
        int f = queue.front();
        queue.pop();
        if (std::find(af.begin(), af.end(), f) != af.end())
          minimum = std::min(minimum, distance[f]);
        if (distance[f] >= std::min(10, minimum + 2))
          continue;
        for (const auto &e : neighbors[f])
          if (distance[e.first] < 0) {
            distance[e.first] = distance[f] + 1;
            queue.push(e.first);
          }
      }
      require(minimum < 100000,
              "Curve samples are separated by a gap at spline " +
                  std::to_string(seg.spline) +
                  "; move "
                  "the curve closer to the surface");
      std::vector<Piece> bestPieces;
      std::vector<int> bestFaces;
      double bestLength = std::numeric_limits<double>::infinity();
      int attempts = 0;
      auto attempt = [&](const std::vector<int> &path) {
        if (path.size() == 2) {
          std::pair<int, int> portal;
          for (const auto &edge : neighbors[path[0]])
            if (edge.first == path[1]) {
              portal = edge.second;
              break;
            }
          Vec3 a = mesh.points[portal.first],
               axis = unit(mesh.points[portal.second] - a);
          // A non-planar quad has no exact common face plane. For
          // this two-face corridor use the planes through each
          // projected sample and the shared edge itself.
          Vec3 np = (p - a).cross(axis), nq = (q - a).cross(axis);
          if (np.norm() > tiny && nq.norm() > tiny) {
            np.normalize();
            nq.normalize();
            if (np.dot(normal(mesh, path[0])) < 0)
              np = -np;
            if (nq.dot(normal(mesh, path[1])) < 0)
              nq = -nq;
            double angle = std::atan2(axis.dot(nq.cross(np)), nq.dot(np));
            Vec3 unfolded = a + Eigen::AngleAxisd(angle, axis) * (q - a);
            Vec3 u = axis, v = np.cross(axis);
            Eigen::Vector2d d((unfolded - p).dot(u), (unfolded - p).dot(v)),
                x((a - p).dot(u), (a - p).dot(v)),
                e((mesh.points[portal.second] - a).dot(u),
                  (mesh.points[portal.second] - a).dot(v));
            double den = cross2(d, e);
            if (std::abs(den) > tiny) {
              double beta = cross2(x, e) / den, t = cross2(x, d) / den;
              if (beta >= -1e-6 && beta <= 1.000001 && t >= -1e-6 &&
                  t <= 1.000001 && (unfolded - p).norm() < bestLength) {
                Vec3 cut = a + std::clamp(t, 0., 1.) *
                                   (mesh.points[portal.second] - a);
                bestLength = (unfolded - p).norm();
                bestPieces.clear();
                bestFaces.clear();
                if ((cut - p).norm() > eps * .01) {
                  bestPieces.push_back(
                      {p, cut, sa.att, attachment(cut, path[0]), s, 0, beta});
                  bestFaces.push_back(path[0]);
                }
                if ((q - cut).norm() > eps * .01) {
                  bestPieces.push_back(
                      {cut, q, attachment(cut, path[1]), sb.att, s, beta, 1});
                  bestFaces.push_back(path[1]);
                }
                return;
              }
            }
          }
        }
        std::vector<Mat3> rotations(path.size(), Mat3::Identity());
        std::vector<Vec3> translations(path.size(), Vec3::Zero());
        std::vector<std::pair<int, int>> portals;
        for (size_t k = 1; k < path.size(); ++k) {
          std::pair<int, int> portal;
          for (const auto &e : neighbors[path[k - 1]])
            if (e.first == path[k]) {
              portal = e.second;
              break;
            }
          portals.push_back(portal);
          Vec3 a = mesh.points[portal.first],
               axis = unit(mesh.points[portal.second] - a),
               na = normal(mesh, path[k - 1]), nb = normal(mesh, path[k]);
          double angle = std::atan2(axis.dot(nb.cross(na)), nb.dot(na));
          Mat3 rotation = Eigen::AngleAxisd(angle, axis).toRotationMatrix();
          rotations[k] = rotations[k - 1] * rotation;
          translations[k] =
              translations[k - 1] + rotations[k - 1] * (a - rotation * a);
        }
        Vec3 target = rotations.back() * q + translations.back(),
             n = normal(mesh, path[0]);
        Vec3 delta = target - p;
        if (delta.norm() < tiny || delta.norm() >= bestLength)
          return;
        Vec3 u = delta.normalized(), v = n.cross(u);
        auto xy = [&](const Vec3 &x) {
          return Eigen::Vector2d((x - p).dot(u), (x - p).dot(v));
        };
        auto d = xy(target);
        Vec3 current = p;
        double alpha = 0;
        std::vector<Piece> traced;
        std::vector<int> tracedFaces;
        for (size_t k = 0; k < path.size(); ++k) {
          Vec3 next = q;
          double beta = 1;
          if (k < portals.size()) {
            Vec3 a = rotations[k] * mesh.points[portals[k].first] +
                     translations[k],
                 b = rotations[k] * mesh.points[portals[k].second] +
                     translations[k];
            Eigen::Vector2d x = xy(a), e = xy(b) - x;
            double den = cross2(d, e);
            if (std::abs(den) < tiny)
              return;
            beta = cross2(x, e) / den;
            double t = cross2(x, d) / den;
            if (beta < alpha - 1e-5 || beta > 1.00001 || t < -.0001 ||
                t > 1.0001)
              return;
            t = std::clamp(t, 0., 1.);
            next = (1 - t) * mesh.points[portals[k].first] +
                   t * mesh.points[portals[k].second];
          }
          if ((next - current).norm() > eps * .01) {
            traced.push_back({current, next, attachment(current, path[k]),
                              attachment(next, path[k]), s, alpha, beta});
            tracedFaces.push_back(path[k]);
          }
          current = next;
          alpha = beta;
        }
        bestLength = delta.norm();
        bestPieces = std::move(traced);
        bestFaces = std::move(tracedFaces);
      };
      // Near a mesh vertex several equally short dual paths exist. A
      // single arbitrary BFS path is insufficient: test unfolded paths
      // and keep the shortest valid surface corridor instead.
      std::vector<int> path;
      std::function<void(int, int)> walk = [&](int f, int remaining) {
        if (++attempts > 2048)
          return;
        path.push_back(f);
        if (std::find(bf.begin(), bf.end(), f) != bf.end())
          attempt(path);
        else if (remaining > 0)
          for (const auto &e : neighbors[f])
            if (distance[e.first] >= 0 && distance[e.first] <= remaining - 1 &&
                std::find(path.begin(), path.end(), e.first) == path.end())
              walk(e.first, remaining - 1);
        path.pop_back();
      };
      for (int f : af)
        walk(f, std::min(10, minimum + 2));
      if (bestPieces.empty()) {
        // Non-planar polygons can make the face-plane unfolding ambiguous.
        // For directly adjacent faces, minimize the two-segment path through
        // their actual shared edge. Bound the detour to reject projection
        // jumps across folds instead of silently routing around the surface.
        for (int a : af)
          for (const auto &edge : neighbors[a]) {
            int b = edge.first;
            if (std::find(bf.begin(), bf.end(), b) == bf.end())
              continue;
            Vec3 x = mesh.points[edge.second.first],
                 d = mesh.points[edge.second.second] - x;
            auto length = [&](double t) {
              Vec3 cut = x + t * d;
              return (p - cut).norm() + (q - cut).norm();
            };
            double lo = 0, hi = 1;
            for (int j = 0; j < 48; ++j) {
              double u = (2 * lo + hi) / 3, v = (lo + 2 * hi) / 3;
              if (length(u) < length(v))
                hi = v;
              else
                lo = u;
            }
            double t = (lo + hi) * .5, l = length(t);
            if (l > 1.5 * (q - p).norm() + 4 * eps || l >= bestLength)
              continue;
            Vec3 mid = x + t * d;
            double beta = (mid - p).norm() / std::max(tiny, l);
            bestLength = l;
            bestPieces = {{p, mid, sa.att, attachment(mid, a), s, 0, beta},
                          {mid, q, attachment(mid, b), sb.att, s, beta, 1}};
            bestFaces = {a, b};
          }
      }
      if (bestPieces.empty()) {
        // Snapped samples on two edges meeting at a mesh vertex may
        // have no straight unfolded corridor (a cone singularity).
        // Preserve the exact edge route through that vertex instead.
        for (int vertex : mesh.faces[sa.att.face]) {
          if (std::find(mesh.faces[sb.att.face].begin(),
                        mesh.faces[sb.att.face].end(),
                        vertex) == mesh.faces[sb.att.face].end())
            continue;
          Vec3 mid = mesh.points[vertex];
          double length = (mid - p).norm() + (q - mid).norm();
          if (length > 2 * (q - p).norm() + eps || length >= bestLength)
            continue;
          double beta = (mid - p).norm() / std::max(tiny, length);
          bestLength = length;
          bestPieces = {
              {p, mid, sa.att, attachment(mid, sa.att.face), s, 0, beta},
              {mid, q, attachment(mid, sb.att.face), sb.att, s, beta, 1}};
          bestFaces = {sa.att.face, sb.att.face};
        }
      }
      require(
          !bestPieces.empty(),
          "No valid unfolded surface corridor at spline " +
              std::to_string(seg.spline) + ", sample " + std::to_string(seg.a) +
              ", faces " + std::to_string(sa.att.face) + "/" +
              std::to_string(sb.att.face) + ", distance " +
              std::to_string(minimum) + ". Add a knot near the sharp feature");
      for (size_t k = 0; k < bestPieces.size(); ++k)
        pieces[bestFaces[k]].push_back(bestPieces[k]);
    }
    // An edge-aligned cut belongs to both incident original polygons,
    // including the pieces generated by the cone-vertex route above.
    auto originalPieces = pieces;
    for (int f = 0; f < int(originalPieces.size()); ++f)
      for (const auto &piece : originalPieces[f]) {
        const auto &ids = mesh.faces[f];
        for (size_t j = 0; j < ids.size(); ++j) {
          int a = ids[j], b = ids[(j + 1) % ids.size()];
          Vec3 start = mesh.points[a], d = mesh.points[b] - start;
          auto on = [&](const Vec3 &x) {
            double t = (x - start).dot(d) / d.squaredNorm();
            return t >= -1e-7 && t <= 1.0000001 &&
                   (x - (start + std::clamp(t, 0., 1.) * d)).norm() < eps * .1;
          };
          if (on(piece.a) && on(piece.b)) {
            for (int other : edgeFaces.at(key(a, b)))
              if (other != f)
                pieces[other].push_back(piece);
            break;
          }
        }
      }
    // All face-local arrangements share original vertex unknowns and exact
    // sample constraints; halfedge walks also represent open-ended cracks.
    std::vector<Label> vertexLabels(mesh.points.size());
    for (const auto &ps : pieces)
      for (const auto &p : ps)
        for (int side = 0; side < 2; ++side) {
          Vec3 x = side ? p.b : p.a;
          auto a = side ? p.bb : p.aa;
          for (int id : a.vertices)
            if ((mesh.points[id] - x).norm() < eps)
              vertexLabels[id] = {p.seg, side ? p.t1 : p.t0, 0};
        }
    auto constraint = [&](const Label &l) {
      const auto &seg = segments[l.seg];
      return std::vector<std::pair<int, double>>{
          {seg.a * 2 + l.side, 1 - l.alpha}, {seg.b * 2 + l.side, l.alpha}};
    };
    for (int f = 0; f < int(mesh.faces.size()); ++f) {
      Vec3 n = normal(mesh, f), u = n.unitOrthogonal(), v = n.cross(u);
      auto xy = [&](const Vec3 &p) {
        return Eigen::Vector2d(p.dot(u), p.dot(v));
      };
      std::vector<CutNode> nodes;
      std::vector<Edge> edges;
      auto addNode = [&](const Vec3 &p, const Attachment &a, int vertex,
                         Label label) {
        for (int i = 0; i < int(nodes.size()); ++i)
          if ((nodes[i].p - p).norm() < eps) {
            if (label.seg >= 0)
              nodes[i].fallback = label;
            return i;
          }
        nodes.push_back({p, a, vertex, label});
        return int(nodes.size()) - 1;
      };
      const auto &fp = mesh.faces[f];
      for (int id : fp) {
        Attachment a{{id, id, id}, Vec3(1, 0, 0), f, mesh.points[id]};
        addNode(mesh.points[id], a, id, vertexLabels[id]);
      }
      for (const auto &p : pieces[f]) {
        addNode(p.a, p.aa, -1, {p.seg, p.t0, 0});
        addNode(p.b, p.bb, -1, {p.seg, p.t1, 0});
      }
      auto addEdge = [&](int a, int b, int s, double t0, double t1,
                         bool boundary) {
        if (a == b)
          return;
        for (auto &e : edges)
          if (key(e.a, e.b) == key(a, b)) {
            if (s >= 0) {
              if (e.seg >= 0 && e.seg != s) {
                bool adjacent = segments[e.seg].a == segments[s].b ||
                                segments[e.seg].b == segments[s].a;
                bool same = segments[e.seg].spline == segments[s].spline;
                if (adjacent && same &&
                    (nodes[a].p - nodes[b].p).norm() < eps * 4)
                  return;
                throw std::runtime_error(
                    "Overlapping projected curves at splines " +
                    std::to_string(segments[e.seg].spline) + "/" +
                    std::to_string(segments[s].spline) + ", overlap " +
                    std::to_string((nodes[a].p - nodes[b].p).norm()) +
                    ". Refit or split the profiles");
              }
              e.seg = s;
              e.t0 = e.a == a ? t0 : t1;
              e.t1 = e.a == a ? t1 : t0;
            }
            e.boundary = e.boundary || boundary;
            return;
          }
        edges.push_back({a, b, s, t0, t1, boundary});
      };
      auto splitLine = [&](const Vec3 &a, const Vec3 &b, int s, double t0,
                           double t1, bool boundary) {
        Vec3 d = b - a;
        double len2 = d.squaredNorm();
        if (len2 < tiny)
          return;
        std::vector<std::pair<double, int>> ids;
        for (int i = 0; i < int(nodes.size()); ++i) {
          double t = (nodes[i].p - a).dot(d) / len2;
          if (t >= -1e-7 && t <= 1.0000001 &&
              (nodes[i].p - (a + std::clamp(t, 0., 1.) * d)).norm() < eps * 2)
            ids.push_back({std::clamp(t, 0., 1.), i});
        }
        std::sort(ids.begin(), ids.end());
        for (size_t i = 1; i < ids.size(); ++i)
          addEdge(ids[i - 1].second, ids[i].second, s,
                  t0 + (t1 - t0) * ids[i - 1].first,
                  t0 + (t1 - t0) * ids[i].first, boundary);
      };
      for (size_t i = 0; i < fp.size(); ++i)
        splitLine(mesh.points[fp[i]], mesh.points[fp[(i + 1) % fp.size()]], -1,
                  0, 1, true);
      for (const auto &p : pieces[f])
        splitLine(p.a, p.b, p.seg, p.t0, p.t1, false);
      // Crossing curves must have an explicit common knot. Rejecting is
      // safer than inventing connectivity and non-manifold constraints.
      for (int i = 0; i < int(edges.size()); ++i)
        if (edges[i].seg >= 0)
          for (int j = 0; j < i; ++j)
            if (edges[j].seg >= 0) {
              const auto &a = edges[i];
              const auto &b = edges[j];
              if (a.a == b.a || a.a == b.b || a.b == b.a || a.b == b.b)
                continue;
              Eigen::Vector2d x = xy(nodes[a.a].p), d = xy(nodes[a.b].p) - x,
                              y = xy(nodes[b.a].p), e = xy(nodes[b.b].p) - y;
              double den = cross2(d, e);
              if (std::abs(den) < tiny)
                continue;
              double s = cross2(y - x, e) / den, t = cross2(y - x, d) / den;
              require(!(s > 1e-6 && s < .999999 && t > 1e-6 && t < .999999),
                      "Curves cross inside a face without a shared knot at "
                      "splines " +
                          std::to_string(segments[a.seg].spline) + " and " +
                          std::to_string(segments[b.seg].spline) +
                          "; split and weld them");
            }
      std::vector<std::vector<int>> adj(nodes.size());
      for (int e = 0; e < int(edges.size()); ++e) {
        adj[edges[e].a].push_back(e * 2);
        adj[edges[e].b].push_back(e * 2 + 1);
      }
      auto from = [&](int h) {
        return h % 2 ? edges[h / 2].b : edges[h / 2].a;
      };
      auto to = [&](int h) { return from(h ^ 1); };
      std::vector<bool> active(nodes.size(), false);
      std::queue<int> q;
      for (int e = 0; e < int(edges.size()); ++e)
        if (edges[e].boundary) {
          active[edges[e].a] = active[edges[e].b] = true;
          q.push(edges[e].a);
          q.push(edges[e].b);
        }
      while (!q.empty()) {
        int a = q.front();
        q.pop();
        for (int h : adj[a])
          if (!active[to(h)]) {
            active[to(h)] = true;
            q.push(to(h));
          }
      }
      for (int i = 0; i < int(nodes.size()); ++i) {
        if (!active[i]) {
          ++stats.removedIslands;
          continue;
        }
        auto &a = adj[i];
        std::sort(a.begin(), a.end(), [&](int x, int y) {
          auto dx = xy(nodes[to(x)].p - nodes[i].p),
               dy = xy(nodes[to(y)].p - nodes[i].p);
          return std::atan2(dx.y(), dx.x()) < std::atan2(dy.y(), dy.x());
        });
      }
      std::vector<bool> seen(edges.size() * 2, false);
      for (int start = 0; start < int(seen.size()); ++start) {
        if (seen[start] || !active[from(start)])
          continue;
        std::vector<int> loop;
        int h = start;
        do {
          require(!seen[h], "Inconsistent cut-face connectivity");
          seen[h] = true;
          loop.push_back(h);
          int b = to(h);
          const auto &a = adj[b];
          auto it = std::find(a.begin(), a.end(), h ^ 1);
          require(it != a.end(), "Missing cut halfedge twin");
          h = a[(int(it - a.begin()) + int(a.size()) - 1) % a.size()];
          require(loop.size() <= seen.size(),
                  "Cut-face traversal did not close");
        } while (h != start);
        Vec3 area = Vec3::Zero();
        for (int d : loop)
          area += nodes[from(d)].p.cross(nodes[to(d)].p);
        if (area.dot(n) <= eps * eps)
          continue;
        Face face;
        std::vector<Vec3> polygon;
        for (int k = 0; k < int(loop.size()); ++k) {
          int dart = loop[k], prev = loop[(k + loop.size() - 1) % loop.size()];
          const auto &node = nodes[from(dart)];
          const auto &e = edges[dart / 2];
          const auto &pe = edges[prev / 2];
          Label label = node.fallback;
          if (e.seg >= 0)
            label = {e.seg, dart % 2 ? e.t1 : e.t0, dart % 2 ? 1 : 0};
          else if (pe.seg >= 0)
            label = {pe.seg, prev % 2 ? pe.t0 : pe.t1, prev % 2 ? 1 : 0};
          Corner c;
          c.vertex = node.vertex;
          c.att = node.att;
          if (label.seg >= 0)
            c.c = constraint(label);
          face.corners.push_back(c);
          polygon.push_back(node.p);
        }
        face.L = polygonLaplacian(polygon);
        faces.push_back(std::move(face));
      }
    }
    // Curvenet constraints take precedence at coincident original vertices.
    for (auto &f : faces)
      for (auto &c : f.corners)
        if (c.vertex >= 0 && vertexLabels[c.vertex].seg >= 0 && c.c.empty())
          c.c = constraint(vertexLabels[c.vertex]);
  }
  void assemble() {
    int nv = int(mesh.points.size()), nc = int(samples.size()) * 2;
    std::vector<bool> constrained(nv, false), used(nv, false),
        anchored(nv, false);
    std::vector<std::vector<std::pair<int, double>>> vertexConstraints(nv);
    for (const auto &f : faces) {
      bool hasConstraint = false;
      for (const auto &c : f.corners) {
        hasConstraint = hasConstraint || !c.c.empty();
        if (c.vertex >= 0) {
          used[c.vertex] = true;
          if (!c.c.empty()) {
            constrained[c.vertex] = true;
            vertexConstraints[c.vertex] = c.c;
          }
        }
      }
      if (hasConstraint)
        for (const auto &c : f.corners)
          for (int v : mesh.faces[c.att.face])
            anchored[v] = true;
    }
    for (auto &f : faces)
      for (auto &c : f.corners)
        if (c.vertex >= 0 && constrained[c.vertex] && c.c.empty())
          c.c = vertexConstraints[c.vertex];
    unknown.assign(nv, -1);
    for (int i = 0; i < nv; ++i) {
      require(used[i], "Mesh has isolated vertices");
      if (!constrained[i]) {
        unknown[i] = int(freeVertices.size());
        freeVertices.push_back(i);
      }
    }
    std::vector<std::vector<int>> adj(nv);
    for (const auto &f : mesh.faces)
      for (size_t i = 0; i < f.size(); ++i) {
        adj[f[i]].push_back(f[(i + 1) % f.size()]);
        adj[f[(i + 1) % f.size()]].push_back(f[i]);
      }
    std::vector<bool> seen(nv, false);
    for (int i = 0; i < nv; ++i)
      if (!seen[i]) {
        std::queue<int> q;
        q.push(i);
        seen[i] = true;
        bool fixed = false;
        while (!q.empty()) {
          int a = q.front();
          q.pop();
          fixed = fixed || constrained[a] || anchored[a];
          for (int b : adj[a])
            if (!seen[b]) {
              seen[b] = true;
              q.push(b);
            }
        }
        if (!fixed) {
          pins.push_back(i);
          ++stats.pinnedComponents;
        }
      }
    std::vector<Triplet> aa, bb;
    for (auto &f : faces) {
      int n = int(f.corners.size());
      for (auto &c : f.corners)
        if (c.vertex >= 0 && c.c.empty())
          c.unknown = unknown[c.vertex];
      for (int i = 0; i < n; ++i)
        if (f.corners[i].unknown >= 0)
          for (int j = 0; j < n; ++j) {
            int row = f.corners[i].unknown;
            double l = f.L(i, j);
            if (f.corners[j].unknown >= 0)
              aa.emplace_back(row, f.corners[j].unknown, l);
            else
              for (auto w : f.corners[j].c)
                bb.emplace_back(row, w.first, -l * w.second);
          }
    }
    for (int v : pins)
      aa.emplace_back(unknown[v], unknown[v], 1);
    int nu = int(freeVertices.size());
    A.resize(nu, nu);
    A.setFromTriplets(aa.begin(), aa.end());
    B.resize(nu, nc);
    B.setFromTriplets(bb.begin(), bb.end());
    if (nu) {
      factor.compute(A);
      require(factor.info() == Eigen::Success,
              "Cut-aware Laplacian factorization failed");
      require((factor.vectorD().array() > 0).all(),
              "Cut-aware Laplacian is not positive definite");
    }
    stats.cutFaces = int(faces.size());
    stats.unknowns = nu;
    for (const auto &f : faces)
      stats.corners += int(f.corners.size());
  }
};

Solver::Solver() : impl_(std::make_unique<Impl>()) {}
Solver::~Solver() = default;
Solver::Solver(Solver &&) noexcept = default;
Solver &Solver::operator=(Solver &&) noexcept = default;
void Solver::bind(const Mesh &mesh, const Net &net, const Options &options) {
  auto start = Clock::now();
  validate(mesh, net);
  require(options.samplesPerEdge >= 1 && std::isfinite(options.samplesPerEdge),
          "Invalid sampling density");
  auto p = std::make_unique<Impl>();
  p->mesh = mesh;
  p->net = net;
  p->options = options;
  BVH bvh(p->mesh);
  Vec3 lo = mesh.points[0], hi = lo;
  double sum = 0;
  int edges = 0;
  for (const auto &v : mesh.points) {
    lo = lo.cwiseMin(v);
    hi = hi.cwiseMax(v);
  }
  p->eps = std::max(1e-10, (hi - lo).norm() * 1e-5);
  for (const auto &f : mesh.faces)
    for (size_t i = 0; i < f.size(); ++i) {
      sum += (mesh.points[f[i]] - mesh.points[f[(i + 1) % f.size()]]).norm();
      ++edges;
    }
  double mean = sum / edges;
  for (const auto &cp : net.points)
    p->cpAttachments.push_back(bvh.query(cp));
  auto snapped = [&](Attachment a) {
    const auto &face = mesh.faces[a.face];
    for (size_t j = 0; j < face.size(); ++j) {
      int id = face[j];
      if ((a.projected - mesh.points[id]).norm() < p->eps) {
        a.vertices = {id, face[(j + 1) % face.size()],
                      face[(j + face.size() - 1) % face.size()]};
        a.bary = Vec3(1, 0, 0);
        a.projected = mesh.points[id];
        return a;
      }
    }
    for (size_t j = 0; j < face.size(); ++j) {
      int ia = face[j], ib = face[(j + 1) % face.size()];
      Vec3 x = mesh.points[ia], d = mesh.points[ib] - x;
      double t = std::clamp((a.projected - x).dot(d) / d.squaredNorm(), 0., 1.);
      Vec3 q = x + t * d;
      if ((q - a.projected).norm() < p->eps) {
        int third = face[(j + 2) % face.size()];
        a.vertices = {ia, ib, third};
        a.bary = Vec3(1 - t, t, 0);
        a.projected = q;
        return a;
      }
    }
    return a;
  };
  for (int s = 0; s < int(net.splines.size()); ++s) {
    std::array<Vec3, 4> cps;
    const auto &sp = net.splines[s];
    double len = 0;
    for (int i = 0; i < 4; ++i) {
      cps[i] = net.points[sp[i]];
      if (i)
        len += (cps[i] - cps[i - 1]).norm();
    }
    int count =
        std::max(2, int(std::ceil(options.samplesPerEdge * len / mean)));
    require(
        int(p->samples.size()) + count + 1 <= options.maxSamples,
        "Too many curve samples; reduce sampling density or split the setup");
    auto ts = arcParameters(cps, count);
    std::vector<int> ids;
    for (double t : ts) {
      Vec3 q = bezier(cps, t);
      ids.push_back(int(p->samples.size()));
      p->samples.push_back({q, snapped(bvh.query(q))});
    }
    for (int i = 0; i < count; ++i)
      p->segments.push_back({ids[i], ids[i + 1], s});
    p->splineSamples.push_back(std::move(ids));
  }
  p->stats.samples = int(p->samples.size());
  p->makeCurves();
  p->cutting(bvh);
  p->assemble();
  p->stats.bindingMs = ms(start);
  impl_ = std::move(p);
}
std::vector<Vec3> Solver::warpedControls(const std::vector<Vec3> &rest) const {
  return impl_->warp(rest);
}
std::vector<Mat3>
Solver::controlRotations(const std::vector<Vec3> &rest) const {
  require(!impl_->faces.empty(), "Solver is not bound");
  require(rest.size() == impl_->mesh.points.size(), "Point count mismatch");
  for (const auto &v : rest)
    require(v.allFinite(), "Non-finite rest surface");
  std::vector<Mat3> out;
  for (const auto &a : impl_->cpAttachments)
    out.push_back(attachmentRotation(impl_->mesh, rest, a));
  return out;
}
const std::vector<Attachment> &Solver::controlAttachments() const {
  return impl_->cpAttachments;
}
const Statistics &Solver::stats() const { return impl_->stats; }
std::vector<Vec3> Solver::solve(const std::vector<Vec3> &rest,
                                const std::vector<Vec3> &posedCP) {
  auto start = Clock::now();
  auto &p = *impl_;
  require(!p.faces.empty(), "Solver is not bound");
  require(rest.size() == p.mesh.points.size() &&
              posedCP.size() == p.net.points.size(),
          "Point count mismatch");
  for (const auto &v : rest)
    require(v.allFinite(), "Non-finite rest surface");
  for (const auto &v : posedCP)
    require(v.allFinite(), "Non-finite posed curve");
  auto base = p.warp(rest);
  std::vector<Vec3> rq(p.samples.size()), pq(rq.size()), rp(rq.size());
  for (int i = 0; i < int(rq.size()); ++i) {
    const auto &a = p.samples[i].att;
    Mat3 F = triangleBasis(rest, a) * triangleBasis(p.mesh.points, a).inverse();
    rp[i] = attached(rest, a);
    rq[i] = rp[i] + F * (p.samples[i].q - a.projected);
  }
  for (int s = 0; s < int(p.net.splines.size()); ++s) {
    std::array<Vec3, 4> b, q;
    for (int j = 0; j < 4; ++j) {
      int v = p.net.splines[s][j];
      b[j] = base[v];
      q[j] = posedCP[v];
    }
    int n = int(p.splineSamples[s].size()) - 1;
    auto bt = arcParameters(b, n), qt = arcParameters(q, n);
    for (int j = 0; j <= n; ++j) {
      int i = p.splineSamples[s][j];
      pq[i] = rq[i] + bezier(q, qt[j]) - bezier(b, bt[j]);
    }
  }
  auto sg = p.gradients(rq, pq);
  std::vector<SideFrame> sampleF(p.samples.size());
  std::vector<std::array<int, 2>> counts(p.samples.size(), {0, 0});
  for (auto &f : sampleF) {
    f.plus.setZero();
    f.minus.setZero();
  }
  for (const auto &c : p.curves)
    for (int j = 0; j < int(c.seg.size()); ++j) {
      int s = c.seg[j];
      const auto &e = p.segments[s];
      for (int i : {e.a, e.b}) {
        sampleF[i].plus += sg[s].plus;
        sampleF[i].minus += sg[s].minus;
        ++counts[i][0];
      }
    }
  int nc = int(p.samples.size()) * 2, nu = int(p.freeVertices.size());
  Eigen::MatrixXd fc(nc, 9), xc(nc, 3);
  for (int i = 0; i < int(p.samples.size()); ++i)
    for (int side = 0; side < 2; ++side) {
      Mat3 f = (side ? sampleF[i].minus : sampleF[i].plus) /
               std::max(1, counts[i][0]);
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
          fc(i * 2 + side, a * 3 + b) = f(a, b);
      xc.row(i * 2 + side) = (pq[i] - f * (rq[i] - rp[i])).transpose();
    }
  Eigen::MatrixXd rhsF = p.B * fc;
  for (int v : p.pins)
    for (int a = 0; a < 3; ++a)
      rhsF(p.unknown[v], a * 3 + a) += 1;
  Eigen::MatrixXd fv(nu, 9);
  if (nu) {
    fv = p.factor.solve(rhsF);
    require(p.factor.info() == Eigen::Success && fv.allFinite(),
            "Gradient solve failed");
  }
  Eigen::MatrixXd rhs = p.B * xc;
  std::vector<Vec3> out(rest.size(), Vec3::Zero());
  std::vector<int> outCount(rest.size(), 0);
  for (const auto &face : p.faces) {
    int n = int(face.corners.size());
    Eigen::MatrixXd Y(n, 3), correction = Eigen::MatrixXd::Zero(n, 3);
    Mat3 F = Mat3::Zero();
    for (const auto &c : face.corners) {
      Eigen::RowVectorXd row = Eigen::RowVectorXd::Zero(9);
      if (c.unknown >= 0)
        row = fv.row(c.unknown);
      else
        for (auto w : c.c)
          row += w.second * fc.row(w.first);
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
          F(a, b) += row[a * 3 + b] / n;
    }
    for (int i = 0; i < n; ++i) {
      const auto &c = face.corners[i];
      Vec3 x = attached(rest, c.att);
      Y.row(i) = (F * x).transpose();
      if (c.unknown < 0) {
        Vec3 interpolated = Vec3::Zero();
        Mat3 cf = Mat3::Zero();
        for (auto w : c.c) {
          interpolated += w.second * rp[w.first / 2];
          for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
              cf(a, b) += w.second * fc(w.first, a * 3 + b);
        }
        correction.row(i) = (cf * (x - interpolated)).transpose();
      }
    }
    // A traced intersection on a folded surface need not lie on the chord
    // joining projected samples. Transport this additional residual too,
    // preserving an exact identity pose instead of introducing bind drift.
    Eigen::MatrixXd ly = face.L * (Y - correction);
    for (int i = 0; i < n; ++i) {
      const auto &c = face.corners[i];
      if (c.unknown >= 0)
        rhs.row(c.unknown) += ly.row(i);
      else if (c.vertex >= 0) {
        Vec3 x = correction.row(i).transpose();
        for (auto w : c.c)
          x += w.second * xc.row(w.first).transpose();
        out[c.vertex] += x;
        ++outCount[c.vertex];
      }
    }
  }
  for (int v : p.pins)
    rhs.row(p.unknown[v]) += rest[v].transpose();
  Eigen::MatrixXd xv(nu, 3);
  if (nu) {
    xv = p.factor.solve(rhs);
    require(p.factor.info() == Eigen::Success && xv.allFinite(),
            "Position solve failed");
    p.stats.relativeResidual =
        (p.A * xv - rhs).norm() / std::max(1.0, rhs.norm());
    require(p.stats.relativeResidual < 1e-6,
            "Position solve is ill-conditioned; inspect degenerate cut cells");
  }
  for (int i = 0; i < int(out.size()); ++i) {
    if (p.unknown[i] >= 0)
      out[i] = xv.row(p.unknown[i]).transpose();
    else {
      require(outCount[i] > 0, "Constrained vertex missing from cut faces");
      out[i] /= outCount[i];
    }
  }
  p.stats.solveMs = ms(start);
  return out;
}
} // namespace cn
