# Created by Rui-Liang Lyu, March 20, 2020
#
# Read out data from a PDB file and organize them in a Pandas DataFrame. Water,
# metal ions, unnatural amino acids are removed in the process.
#

import pandas as pd
from urllib.request import urlopen
from tqdm import tqdm


natural_aas = [
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL"
]

columns = [
    "chainid",
    "resid",
    "resname",
    "atomid",
    "atomname",
    "x",
    "y",
    "z",
    "element",
    "altloc"
]

parse_atom_line = lambda line: {
    "chainid":  line[21],
    "resid":    line[22:27].strip(" "),
    "resname":  line[17:20],
    "atomid":   int(line[6:11]),
    "atomname": line[12:16].strip(" "),
    "x":        float(line[30:38]),
    "y":        float(line[38:46]),
    "z":        float(line[47:54]),
    "element":  line[76:78].strip(" "),
    "altloc":   line[16]
}


def read_local_pdb_file(path):
    with open(path, "r") as f:
        lines = f.readlines()
    return lines


def fetch_pdb_from_rcsb(pdbid):
    url = "https://files.rcsb.org/download/{}.pdb".format(pdbid)
    lines = [line.decode("utf-8") for line in urlopen(url)]
    return lines


def parse(pdb, verbose=False):
    try:
        lines = read_local_pdb_file(pdb)
    except FileNotFoundError:
        lines = fetch_pdb_from_rcsb(pdb)
    except Exception:
        print("Error: failed to open local file `{}` or fetch from RCSB PDB".format(pdb))
        return None

    lines = [line for line in lines if line.startswith("ATOM  ")]

    data = pd.DataFrame(columns=columns)

    if verbose:
        iterator = tqdm(lines, ascii=False, desc="parsing {}".format(pdb))
    else:
        iterator = lines

    for line in iterator:
        atom = parse_atom_line(line)

        # change HSE to HIS
        if atom["resname"] == "HSE":
            atom["resname"] = "HIS"

        # remove those that are not natural amino acids
        if atom["resname"] not in natural_aas:
            continue

        # remove hydrogens, metals, and other elements not common in proteins
        if atom["element"] not in ["C", "N", "O", "S"]:
            continue

        # if a residue has altlocs, only keep the one first seen
        if atom["altloc"] != " ":
            altlocs = set(data[data.resid == atom["resid"]].altloc)
            if len(altlocs) > 0 and atom["altloc"] not in altlocs:
                continue

        data = data.append(atom, ignore_index=True)  # note that pd.DataFrame.append method is not in-place

    return data.drop("altloc", axis=1)  # drop the altloc column since it is useless for calculating SASA
