"""Host-only authoring-data checks; these do not run the Maya API."""

import importlib.util
import math
import unittest
from pathlib import Path

module_path = Path(__file__).resolve().parents[1] / "scripts" / "curvenet_maya.py"
spec = importlib.util.spec_from_file_location("curvenet_maya", module_path)
cn = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cn)


class AuthoringDataTests(unittest.TestCase):
    def test_shared_knots_and_control_roles(self):
        bundle = cn.validate_bundle(cn.demo_bundle())
        net = bundle["network"]
        segments = [[net["points"][i] for i in span] for span in net["splines"]]
        graph = cn.graph_from_segments(segments)
        self.assertEqual(len(graph["splines"]), 3)
        self.assertEqual(sum(p < 0 for p in graph["parents"]), 4)
        self.assertEqual(len({span[0] for span in graph["splines"]}), 1)

    def test_cubic_reconstruction(self):
        controls = [[0.0, 0.0, 0.0], [1.0, 2.0, 0.0], [2.0, 1.0, 0.0], [4.0, 0.0, 0.0]]

        def sample(t):
            s = 1 - t
            return [
                s**3 * controls[0][i]
                + 3 * s * s * t * controls[1][i]
                + 3 * s * t * t * controls[2][i]
                + t**3 * controls[3][i]
                for i in range(3)
            ]

        actual = cn._bezier_from_samples([sample(t) for t in (0, 1 / 3, 2 / 3, 1)])
        self.assertLess(max(math.dist(a, b) for a, b in zip(controls, actual)), 1e-12)

    def test_reject_bad_graph(self):
        with self.assertRaises(ValueError):
            cn.graph_from_segments([], 1e-4)
        with self.assertRaises(ValueError):
            cn.graph_from_segments([[[0, 0, 0]] * 4])
        with self.assertRaises(ValueError):
            cn.graph_from_segments([], -1)


if __name__ == "__main__":
    unittest.main()
