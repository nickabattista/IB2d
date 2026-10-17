from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from generate_case import build_case


class GeometryTests(unittest.TestCase):
    def test_generated_ids_spacing_and_metadata(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            case_dir = Path(temporary_directory)
            metadata = build_case(case_dir, nx=32, ny=32, dt=2.5e-5, steps=4)

            vertex_lines = (case_dir / "plate.vertex").read_text().splitlines()
            target_lines = (case_dir / "plate.target").read_text().splitlines()
            count = int(vertex_lines[0])

            self.assertEqual(count, metadata["geometry"]["marker_count"])
            self.assertEqual(len(vertex_lines) - 1, count)
            self.assertEqual(len(target_lines) - 1, count)
            self.assertEqual(
                [int(line.split()[0]) for line in target_lines[1:]],
                list(range(count)),
            )
            self.assertAlmostEqual(
                metadata["geometry"]["effective_quadrature_chord"],
                count * metadata["geometry"]["marker_weight"],
            )
            loaded = json.loads((case_dir / "case_metadata.json").read_text())
            self.assertEqual(loaded["time"]["steps"], 4)


if __name__ == "__main__":
    unittest.main()
