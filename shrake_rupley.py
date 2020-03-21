# Created by Rui-Liang Lyu, March 20, 2020
#
# Calculate solvent accessible surface area (SASA) of each atom in a protein
# by the Shrake-Rupley algorithm.
#

import numpy as np
from numpy import pi, sin, cos, sqrt
import pandas as pd
from tqdm import tqdm

from pdb_parser import parse
from atomic_radii import get_radius


def get_points_on_sphere(center=(0., 0., 0.), radius=1., n=100):
    pts = []

    inc = pi * (3 - sqrt(5))
    off = 2 / n

    for k in range(n):
        y = k * off - 1 + (off / 2)
        r = sqrt(1 - y * y)
        phi = k * inc
        pts.append([cos(phi) * r, y, sin(phi) * r])

    pts = np.array(pts) * radius + np.array(center)
    return pts


def shrake_rupley(atoms, probe_radius=1.4, n_samples=150, verbose=False):
    n_atoms = len(atoms)
    centers = []
    radii = []
    results = []

    for idx, atom in atoms.iterrows():
        centers.append((atom["x"], atom["y"], atom["z"]))
        radii.append(get_radius(atom["resname"], atom["atomname"], atom["element"]))

    if verbose:
        iterator = tqdm(zip(centers, radii),
                        total=n_atoms,
                        ascii=False,
                        desc="calculating surface area")
    else:
        iterator = zip(centers[:1], radii[:1])

    for center, radius in iterator:
        pts = get_points_on_sphere(center=center,
                                   radius=(radius + probe_radius),
                                   n=n_samples)

        pts = pts.repeat(len(atoms), 0).reshape(n_samples, n_atoms, 3)
        d2 = np.sum((pts - np.array(centers)) ** 2, axis=2)
        r2 = (np.array(radii) + probe_radius) ** 2
        r2 = np.stack([r2] * n_samples)
        n_outsiders = np.sum(np.all(d2 >= (r2 * 0.99), axis=1))  # the 0.99 factor to account for numerical errors in the calculation of d2

        area = 4 * pi * ((radius + probe_radius) ** 2) * n_outsiders / n_samples
        results.append(area)

    return results
