"""Writes ref_data_laliberte_2009.json: density and viscosity of single-electrolyte brines
computed by the thermo package from the Laliberte (2009) coefficients, inside each fitted
range. Run with python3 laliberte_2009_reference.py > ref_data_laliberte_2009.json"""
import csv
import json
import os
import sys

import chemicals
import numpy as np
import thermo
import thermo.electrochem as electrochem

NAMES = ["LiCl", "Li2SO4", "SrCl2", "NaBr", "KBr", "NaHCO3", "Na2CO3", "KHCO3", "K2CO3"]
TSV = os.path.join(os.path.dirname(chemicals.__file__), "Electrolytes", "Laliberte2009.tsv")

with open(TSV) as f:
    rows = {r[2]: r for r in csv.reader(f, delimiter="\t")}

out = {"_source": f"thermo {thermo.__version__}, chemicals {chemicals.__version__}"}
for name in NAMES:
    row = rows[name]
    cas = row[1]
    entry = {}
    # Columns: fitted tmin, tmax (C) and wmax for density, then for viscosity
    for prop, cols, fn in (("density", (9, 10, 11), electrochem.Laliberte_density),
                           ("viscosity", (19, 20, 21), electrochem.Laliberte_viscosity)):
        tmin, tmax, wmax = (float(row[c]) for c in cols)
        temps = np.unique(np.linspace(max(tmin, 0.0), tmax, 4))
        fracs = np.linspace(0.1, 0.9, 4) * wmax
        entry[prop] = [[float(t + 273.15), float(w), float(fn(t + 273.15, [w], [cas]))]
                       for t in temps for w in fracs]
    out[name] = entry
json.dump(out, sys.stdout, indent=1)
sys.stdout.write("\n")
