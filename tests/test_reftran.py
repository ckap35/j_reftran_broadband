import unittest

import numpy as np

from j_reftran_broadband import jreftran_rt


class ThinFilmTests(unittest.TestCase):
    def test_single_interface_matches_fresnel(self):
        r, t, R, T, A = jreftran_rt(500, [np.nan, np.nan], [1, 1.5], 0, 0)
        np.testing.assert_allclose([r, t, R, T, A], [-0.2, 0.8, 0.04, 0.96, 0], atol=1e-14)

    def test_original_absorbing_film_example(self):
        r, t, R, T, A = jreftran_rt(
            500, [np.nan, 200, np.nan], [1, 0.9707 + 1.8562j, 1], 0, 0
        )
        np.testing.assert_allclose([r, t], [-0.4622 - 0.5066j, -0.0047 + 0.0097j], atol=6e-5)
        np.testing.assert_allclose([R, T, A], [0.4702, 1.1593e-4, 0.5296], atol=6e-5)
        self.assertAlmostEqual(R + T + A, 1)

    def test_lossless_stack_conserves_power(self):
        for wavelength in [400, 550, 900]:
            for angle in [0, 30, 65]:
                for polarization in [0, 1]:
                    with self.subTest(wavelength=wavelength, angle=angle, pol=polarization):
                        _, _, R, T, A = jreftran_rt(
                            wavelength, [np.nan, 90, 140, np.nan],
                            [1, 1.45, 2.1, 1.5], angle, polarization
                        )
                        self.assertAlmostEqual(R + T, 1)
                        self.assertAlmostEqual(A, 0)

    def test_polarizations_agree_at_normal_incidence(self):
        args = (600, [np.nan, 100, np.nan], [1, 2 + 0.2j, 1.5], 0)
        np.testing.assert_allclose(jreftran_rt(*args, 0), jreftran_rt(*args, 1))

    def test_boundary_thicknesses_are_ignored(self):
        args = ([1, 1.5, 1], 25, 0)
        expected = jreftran_rt(500, [np.nan, 100, np.nan], *args)
        np.testing.assert_allclose(jreftran_rt(500, [0, 100, 999], *args), expected)


if __name__ == "__main__":
    unittest.main()
