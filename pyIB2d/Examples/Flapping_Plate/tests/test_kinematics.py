from __future__ import annotations

import json
import math
import os
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from case_config import (
    MU,
    NOMINAL_CHORD,
    PEAK_TRANSLATION_SPEED,
    PERIOD,
    REYNOLDS_NUMBER,
    RHO,
    STROKE_AMPLITUDE,
    default_motion_metadata,
    marker_spacing,
    reynolds_number,
    viscosity_for_reynolds,
)
from kinematics import (
    FlappingMotion,
    body_coordinates_from_world,
    peak_translation_speed,
    point_velocities,
    pose_at_time,
    velocity_at_time,
    world_coordinates_from_body,
)
from update_Target_Point_Positions import (
    reset_reference,
    update_Target_Point_Positions,
)


class ScalingTests(unittest.TestCase):
    def test_fixed_reynolds_scaling(self):
        self.assertAlmostEqual(STROKE_AMPLITUDE / NOMINAL_CHORD, 2.5)
        self.assertAlmostEqual(
            peak_translation_speed(STROKE_AMPLITUDE, PERIOD),
            PEAK_TRANSLATION_SPEED,
        )
        self.assertAlmostEqual(MU, 0.375)
        self.assertAlmostEqual(
            reynolds_number(RHO, PEAK_TRANSLATION_SPEED, NOMINAL_CHORD, MU),
            REYNOLDS_NUMBER,
        )
        self.assertAlmostEqual(
            viscosity_for_reynolds(
                RHO, PEAK_TRANSLATION_SPEED, NOMINAL_CHORD, REYNOLDS_NUMBER
            ),
            MU,
        )

    def test_marker_spacing_is_half_the_fluid_spacing(self):
        self.assertAlmostEqual(marker_spacing(1.6, 512), 0.5 * 1.6 / 512)

    def test_peak_speed_rejects_nonpositive_period(self):
        with self.assertRaises(ValueError):
            peak_translation_speed(STROKE_AMPLITUDE, 0.0)


class PoseTests(unittest.TestCase):
    def setUp(self):
        self.motion = FlappingMotion.from_mapping(default_motion_metadata())
        self.beta = math.radians(self.motion.stroke_angle_degrees)
        self.alpha_mean = math.radians(
            self.motion.mean_angle_of_attack_degrees
        )
        self.pitch_amplitude = math.radians(
            self.motion.pitch_amplitude_degrees
        )

    def test_cardinal_phase_poses_follow_angle_convention(self):
        pose0 = pose_at_time(0.0, self.motion)
        np.testing.assert_allclose(
            pose0.center,
            np.asarray(self.motion.stroke_midpoint)
            + 0.5
            * self.motion.stroke_amplitude
            * np.asarray((math.cos(self.beta), math.sin(self.beta))),
            atol=1.0e-15,
        )
        self.assertAlmostEqual(pose0.alpha_radians, self.alpha_mean)
        self.assertAlmostEqual(
            pose0.theta_radians, self.beta - self.alpha_mean
        )

        quarter = pose_at_time(0.25 * self.motion.period, self.motion)
        np.testing.assert_allclose(
            quarter.center, self.motion.stroke_midpoint, atol=1.0e-15
        )
        self.assertAlmostEqual(
            quarter.alpha_radians, self.alpha_mean + self.pitch_amplitude
        )
        self.assertAlmostEqual(
            quarter.theta_radians,
            self.beta - self.alpha_mean - self.pitch_amplitude,
        )

        half = pose_at_time(0.5 * self.motion.period, self.motion)
        np.testing.assert_allclose(
            half.center,
            np.asarray(self.motion.stroke_midpoint)
            - 0.5
            * self.motion.stroke_amplitude
            * np.asarray((math.cos(self.beta), math.sin(self.beta))),
            atol=1.0e-15,
        )
        self.assertAlmostEqual(half.alpha_radians, self.alpha_mean)

    def test_velocity_matches_centered_stroke_and_theta_rate(self):
        quarter = velocity_at_time(0.25 * self.motion.period, self.motion)
        expected = -peak_translation_speed(
            self.motion.stroke_amplitude, self.motion.period
        ) * np.asarray((math.cos(self.beta), math.sin(self.beta)))
        np.testing.assert_allclose(quarter.center_velocity, expected, atol=1.0e-13)
        self.assertAlmostEqual(quarter.theta_rate, 0.0, places=12)

        zero = velocity_at_time(0.0, self.motion)
        self.assertAlmostEqual(zero.center_velocity_x, 0.0)
        self.assertAlmostEqual(zero.center_velocity_y, 0.0)
        self.assertAlmostEqual(zero.theta_rate, -zero.alpha_rate)

    def test_midchord_transform_is_rigid_and_invertible(self):
        body = np.asarray(((-0.05, 0.0), (0.0, 0.0), (0.05, 0.0)))
        world = world_coordinates_from_body(body, 0.37 * self.motion.period, self.motion)
        recovered = body_coordinates_from_world(
            world, 0.37 * self.motion.period, self.motion
        )
        np.testing.assert_allclose(recovered, body, atol=1.0e-15)
        np.testing.assert_allclose(
            world[1],
            pose_at_time(0.37 * self.motion.period, self.motion).center,
            atol=1.0e-15,
        )
        self.assertAlmostEqual(np.linalg.norm(world[2] - world[0]), 0.1)

    def test_analytic_point_velocity_matches_finite_difference(self):
        body = np.asarray(((-0.05, 0.0), (0.03, 0.0)))
        time = 0.173 * self.motion.period
        # Large enough to avoid subtractive cancellation while remaining far
        # below any resolved fraction of the sinusoidal period.
        epsilon = 1.0e-6 * self.motion.period
        finite_difference = (
            world_coordinates_from_body(body, time + epsilon, self.motion)
            - world_coordinates_from_body(body, time - epsilon, self.motion)
        ) / (2.0 * epsilon)
        np.testing.assert_allclose(
            point_velocities(body, time, self.motion),
            finite_difference,
            rtol=2.0e-8,
            atol=2.0e-8,
        )


class TargetUpdateTests(unittest.TestCase):
    def test_midpoint_update_is_absolute_and_drift_free(self):
        motion_metadata = default_motion_metadata()
        motion = FlappingMotion.from_mapping(motion_metadata)
        body = np.column_stack((np.linspace(-0.05, 0.05, 9), np.zeros(9)))
        initial = world_coordinates_from_body(body, 0.0, motion)
        targets = np.column_stack(
            (np.arange(body.shape[0]), initial, np.full(body.shape[0], 1.0e9))
        )

        with tempfile.TemporaryDirectory() as temporary_directory:
            case_dir = Path(temporary_directory)
            (case_dir / "case_metadata.json").write_text(
                json.dumps({"motion": motion_metadata}), encoding="utf-8"
            )
            previous_directory = Path.cwd()
            try:
                os.chdir(case_dir)
                reset_reference()
                dt = 1.25e-6
                time = 0.31 * motion.period
                first = update_Target_Point_Positions(dt, time, targets.copy())
                repeated = update_Target_Point_Positions(dt, time, first.copy())
            finally:
                os.chdir(previous_directory)
                reset_reference()

        expected = world_coordinates_from_body(body, time + 0.5 * dt, motion)
        np.testing.assert_allclose(first[:, 1:3], expected, atol=1.0e-15)
        np.testing.assert_allclose(repeated, first, atol=1.0e-15)


if __name__ == "__main__":
    unittest.main()
