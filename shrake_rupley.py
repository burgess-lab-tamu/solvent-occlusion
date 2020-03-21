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


def calculate_surface_area(atoms, probe_radius=1.4, n_samples=150, verbose=False):
    n_atoms = len(atoms)
    centers = []
    radii = []
    areas = []

    for idx, atom in atoms.iterrows():
        centers.append((atom["x"], atom["y"], atom["z"]))
        radii.append(get_radius(atom["resname"], atom["atomname"], atom["element"]))

    if verbose:
        iterator = tqdm(zip(centers, radii),
                        total=n_atoms,
                        ascii=False,
                        desc="calculating surface area")
    else:
        iterator = zip(centers, radii)

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
        areas.append(area)

    return pd.concat([
        atoms.chainid,
        atoms.resid,
        atoms.resname,
        pd.Series(areas, index=atoms.index, name="area")
    ], axis=1)


def parse_results_by_residue(results):
    resinfo = []
    areas = []

    for idx, row in results.iterrows():
        res = (row["chainid"], row["resid"], row["resname"])
        if res not in resinfo:
            resinfo.append(res)

    for res in resinfo:
        area = results[(results.chainid == res[0]) &
                       (results.resid == res[1]) &
                       (results.resname == res[2])].area
        areas.append(sum(area))

    resinfo = np.array(resinfo)

    return pd.DataFrame({
        "chainid": resinfo[:, 0],
        "resid": resinfo[:, 1],
        "resname": resinfo[:, 2],
        "area": areas
    })
