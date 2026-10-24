"""Reflection and transmission of nonmagnetic thin-film stacks.

Adapted from Shawn Divitt's MATLAB jreftran calculator:
https://www.mathworks.com/matlabcentral/fileexchange/50923
"""

import argparse
from pathlib import Path

import numpy as np


def sind(angle):
    """Sine in degrees, with exact zeros at scalar multiples of 180."""
    value = np.sin(np.radians(angle))
    turns = angle / 180.0
    if np.ndim(turns) == 0 and np.isfinite(turns) and turns == np.trunc(turns):
        return 0.0
    return value


def jreftran_rt(wavelength, d, n, t0, polarization):
    """Return amplitude coefficients r, t and power fractions R, T, A.

    Wavelength and thicknesses d are in nm; t0 is in degrees. Layers run from
    incident medium to substrate. The first and last thicknesses are ignored
    and may be np.nan. Refractive indices use n + i*k for absorption.
    Use polarization=0 for TE (s) and polarization=1 for TM (p).
    Power fractions assume lossless incident and exit media with propagating
    waves.
    """
    d = np.asarray(d, dtype=float)
    n = np.asarray(n, dtype=complex)
    if d.ndim != 1 or n.ndim != 1 or d.shape != n.shape or len(d) < 2:
        raise ValueError("d and n must be one-dimensional arrays of equal length >= 2")

    impedance = 376.730313  # Free-space impedance, ohms.
    admittance = n / impedance
    sin_theta = n[0] / n * sind(t0)
    cos_theta = np.sqrt(1 - sin_theta**2)
    if polarization == 0:
        eta = admittance * cos_theta
    else:
        eta = admittance / cos_theta

    matrix = np.identity(2, dtype=complex)
    for j in range(1, len(d) - 1):
        # Characteristic matrix of the finite film.
        phase = -2 * np.pi * n[j] * d[j] * cos_theta[j] / wavelength
        layer = np.array([
            [np.cos(phase), 1j * np.sin(phase) / eta[j]],
            [1j * eta[j] * np.sin(phase), np.cos(phase)],
        ])
        matrix = matrix @ layer

    incident, exit_medium = eta[0], eta[-1]
    field = matrix[0, 0] + matrix[0, 1] * exit_medium
    flux = matrix[1, 0] + matrix[1, 1] * exit_medium
    denominator = incident * field + flux

    r = (incident * field - flux) / denominator
    t = 2 * incident / denominator
    R = abs(r)**2
    T = ((exit_medium.real / incident) * abs(t)**2).real
    A = (4 * incident * (field * np.conj(flux) - exit_medium).real
         / abs(denominator)**2).real
    return r, t, R, T, A


def main(argv=None):
    parser = argparse.ArgumentParser(description="Plot a 20-film wavelength sweep.")
    parser.add_argument("--output", type=Path, help="Save the figure instead of opening a window")
    args = parser.parse_args(argv)

    import matplotlib.pyplot as plt

    # Twenty absorbing films, each preceded by a 2.4 um air gap.
    d = np.r_[np.nan, np.tile([2400.0, 10.0], 20), np.nan]
    n = np.r_[1.0, np.tile([1.0, 0.9707 + 1.8562j], 20), 1.0]
    wavelengths = np.logspace(2, 4, num=1000)
    fractions = np.array([
        jreftran_rt(wavelength, d, n, 0, 0)[2:] for wavelength in wavelengths
    ])

    fig, ax = plt.subplots(figsize=(8, 4.5))
    for values, label in zip(fractions.T, ["Reflection", "Transmission", "Absorption"]):
        ax.semilogx(wavelengths, values, label=label)
    ax.set_xlabel("Wavelength [nm]")
    ax.set_ylabel("Power fraction")
    ax.set_title("20-film stack with fixed refractive indices")
    ax.legend()
    fig.tight_layout()
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.output, dpi=160)
        plt.close(fig)
    else:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
