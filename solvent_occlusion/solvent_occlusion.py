# Created by Rui-Liang Lyu, March 21, 2020
#

import multiprocessing as mp
import pandas as pd

from .shrake_rupley import shrake_rupley, shrake_rupley_mp


def get_areas(atoms, chainids):
    areas_all = shrake_rupley(atoms, by_residue=True)

    areas_sgl = []
    for chainid in chainids:
        chain = atoms[atoms.chainid == chainid]
        areas_sgl.append(shrake_rupley(chain, by_residue=True))

    areas_sgl = pd.concat(areas_sgl).reset_index()
    return areas_all, areas_sgl


def get_areas_mp(atoms, chainids):
    queue = mp.Queue()
    processes = [mp.Process(target=shrake_rupley_mp,
                            args=(atoms, queue, "all"),
                            kwargs={"by_residue": True})]

    for chainid in chainids:
        chain = atoms[atoms.chainid == chainid]
        p = mp.Process(target=shrake_rupley_mp,
                       args=(chain, queue, chainid),
                       kwargs={"by_residue": True})
        processes.append(p)

    for p in processes:
        p.start()

    for p in processes:
        p.join()

    outputs = [queue.get() for p in processes]
    outputs = {output[0]: output[1] for output in outputs}

    areas_all = outputs["all"]
    areas_sgl = [outputs[chainid] for chainid in chainids]
    areas_sgl = pd.concat(areas_sgl).reset_index()
    return areas_all, areas_sgl


def sort_areas(areas_all, areas_sgl):
    indices = [0] * len(areas_sgl)
    for idx, row in areas_sgl.iterrows():
        indices[idx] = areas_all.loc[(areas_all.chainid == row["chainid"]) &
                                     (areas_all.resid == row["resid"])].index[0]

    areas_sgl.index = indices
    areas_sgl = areas_sgl.sort_index()

    return areas_sgl


def get_occlusions(atoms, chainids=None, use_mp=False):
    if chainids is None:
        chainids = sorted(list(set(atoms.chainid)))
    atoms = atoms[atoms.chainid.isin(chainids)]

    if use_mp:
        areas_all, areas_sgl = get_areas_mp(atoms, chainids)
    else:
        areas_all, areas_sgl = get_areas(atoms, chainids)

    areas_sgl = sort_areas(areas_all, areas_sgl)

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
