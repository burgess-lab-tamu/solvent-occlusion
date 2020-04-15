# Created by Rui-Liang Lyu, March 21, 2020
#

import pandas as pd

from .shrake_rupley import get_areas


def get_occlusions(atoms, chainids=None, verbose=False):
    if chainids is None:
        chainids = sorted(list(set(atoms.chainid)))

    atoms = atoms[atoms.chainid.isin(chainids)]
    areas_all = get_areas(atoms,
                          by_residue=True,
                          verbose=verbose,
                          _desc="Calc areas of PPI complex")

    areas_sgl = []
    for chainid in chainids:
        chain = atoms[atoms.chainid == chainid]
        areas_sgl.append(get_areas(chain,
                                   by_residue=True,
                                   verbose=verbose,
                                   _desc="Calc areas of chain {}".format(chainid)))

    areas_sgl = pd.concat(areas_sgl).reset_index()

    # Important: sort areas_sgl so that it is aligned with areas_all
    indices = [0] * len(areas_sgl)
    for idx, row in areas_sgl.iterrows():
        indices[idx] = areas_all.loc[(areas_all.chainid == row["chainid"]) &
                                     (areas_all.resid == row["resid"])].index[0]

    areas_sgl.index = indices
    areas_sgl = areas_sgl.sort_index()

    occlusions = []
    for area_sgl, area_all in zip(areas_sgl.area, areas_all.area):
        if area_sgl == 0:
            occlusions.append(0)
        else:
            occlusions.append((area_sgl - area_all) / area_sgl)

    return pd.DataFrame({
        "chainid": areas_sgl.chainid.tolist(),
        "resid": areas_sgl.resid.tolist(),
        "resname": areas_sgl.resname.tolist(),
        "area_sgl": areas_sgl.area.tolist(),
        "area_all": areas_all.area.tolist(),
        "occlusion": occlusions
    })
