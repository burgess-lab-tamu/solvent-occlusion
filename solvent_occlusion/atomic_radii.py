# Created by Rui-Liang Lyu, March 20, 2020
#
# A function that reads the resname and atomname of an atoms, then returns its
# atomic radius suggested by Shrake & Rupley, 1973.
#

radii = {
    "C":      2.0,   # aliphatic carbons; default if unknown
    "C_aro":  1.85,  # aromatic carbons
    "C_carb": 1.5,   # carbons in carbonyls
    "N":      1.5,
    "O":      1.4,
    "S":      1.85
}

arg = {
    "CG": "C",
    "CD": "C",
    "CZ": "C_carb"
}

asn = {"CG": "C_carb"}

gln = {
    "CG": "C",
    "CD": "C_carb"
}

his = {
    "CE1": "C_aro",
    "AE1": "C_aro",
    "CD2": "C_aro",
    "AD2": "C_aro",
    "CG":  "C_aro"
}

ile = {
    "CG1": "C",
    "CG2": "C",
    "CD1": "C"
}

leu = {
    "CG":  "C",
    "CD1": "C",
    "CD2": "C"
}

lys = {
    "CG": "C",
    "CD": "C",
    "CE": "C"
}

met = {
    "CG": "C",
    "CE": "C"
}

phe = {
    "CG":  "C_aro",
    "CD1": "C_aro",
    "CD2": "C_aro",
    "CE1": "C_aro",
    "CE2": "C_aro",
    "CZ":  "C_aro"
}

pro = {
    "CG": "C",
    "CD": "C"
}

thr = {
    "CB":  "C",
    "CG2": "C"
}

trp = {
    "CG":  "C_aro",
    "CD1": "C_aro",
    "CE2": "C_aro",
    "CZ2": "C_aro",
    "CH2": "C_aro",
    "CZ3": "C_aro",
    "CE3": "C_aro",
    "CD2": "C_aro"
}

val = {
    "CG1": "C",
    "CG2": "C"
}

residues = {
    "ARG": arg, "ASN": asn, "ASP": asn, "GLN": gln, "GLU": gln,
    "HIS": his, "ILE": ile, "LEU": leu, "LYS": lys, "MET": met,
    "PHE": phe, "PRO": pro, "THR": thr, "TRP": trp, "TYR": phe,
    "VAL": val
}


def get_carbon_radius(resname, atomname):
    if atomname in ["CA", "CB"]:
        return radii["C"]
    elif atomname == "C":
        return radii["C_carb"]
    else:
        return radii[residues[resname][atomname]]


def get_radius(resname, atomname, element):
    if element == "C":
        return get_carbon_radius(resname, atomname)
    else:
        return radii[element]
