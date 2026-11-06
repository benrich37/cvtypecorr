"""Step 2: free energy of each elementary step on one surface, before and after the corrections.

The steps are the transitions among *, *NO3, *NO3H, *O, *OH, *NO2 and *H in the step map of the NO3RR pathway
analysis (research: NO3_ads_gen/anl/mnc/no3rr_pathway_map.py), on FeN4-i at 0 and -1 V vs SHE, pH 6.3.

How a step energy is built (see energies.step_energy):
    dG_raw     = sum G(products) - sum G(reactants)
                 G(*X)  = lowest grand free energy of the cached *X configurations at the bias ("*" = clean slab)
                 G(mol) = grand free energy of the free molecule at the bias, + kT ln c for H3O+ (10^-pH M) and
                          H2O (55.34 M)
    correction = sum c(label) over products - over reactants, with each adsorbate labelled by the free species it
                 is taken to share its DFT error with (energies.ADSORBATE_LABEL: *NO3 -> NO3-, *OH -> OH-, ...)
    dG         = dG_raw + correction
A proton is H3O+ -> H2O; its electron is implicit (all G are grand free energies at the electrode potential), so a
"H3O+" step is a proton-coupled electron transfer.

Run: python surface_steps.py   -> output/FeN4-i_steps.csv, output/FeN4-i_steps.md; checks against the research code
"""
import csv

from energies import HERE, REFERENCE, SETTINGS, bias_label, side_str, step_energy
from fit_corrections import fit_corrections

SURF = REFERENCE["steps"]["surf"]
BIASES = ["0V", "n1p0V"]
CONC = SETTINGS["conc_M"]["steps"]

# (step map key, reactants, products)
STEPS = [
    ("* -> NO3", [("*", 1), ("NO3-", 1)], [("*NO3", 1)]),
    ("* -> O", [("*", 1), ("NO3-", 1)], [("*O", 1), ("NO2-", 1)]),
    ("* -> NO2", [("*", 1), ("NO2-", 1)], [("*NO2", 1)]),
    ("* -> H", [("*", 1), ("H3O+", 1)], [("*H", 1), ("H2O", 1)]),
    ("NO3 -> NO3H", [("*NO3", 1), ("H3O+", 1)], [("*NO3H", 1), ("H2O", 1)]),
    ("NO3 -> O", [("*NO3", 1)], [("*O", 1), ("NO2-", 1)]),
    ("NO3H -> NO2", [("*NO3H", 1), ("H3O+", 1)], [("*NO2", 1), ("H2O", 2)]),
    ("O -> OH", [("*O", 1), ("H3O+", 1)], [("*OH", 1), ("H2O", 1)]),
    ("OH -> *", [("*OH", 1), ("H3O+", 1)], [("*", 1), ("H2O", 2)]),
    ("H -> *", [("*H", 1), ("H3O+", 1)], [("*", 1), ("H2O", 1), ("H2", 1)]),
]


def main():
    corrections = fit_corrections()[0].corrections
    research = {(r["step"], r["bias"]): r for r in REFERENCE["steps"]["rows"]}
    rows, max_err = [], 0.0
    for key, reactants, products in STEPS:
        for bias in BIASES:
            s = step_energy(reactants, products, SURF, bias, corrections, CONC)
            ref = research[(key, bias)]
            max_err = max(max_err, abs(s["dG_raw"] - ref["dG_raw"]), abs(s["dG"] - ref["dG_corrected"]))
            rows.append({"step": key, "equation": f"{side_str(reactants)} -> {side_str(products)}", "bias": bias,
                         "U_vs_SHE_V": round(SETTINGS["bias_vs_SHE_V"][bias], 6),
                         "dG_raw_eV": s["dG_raw"], "correction_eV": s["correction"], "dG_eV": s["dG"],
                         "correction_terms": s["correction_terms"] or "0 (labels cancel)",
                         "calc_dirs": "; ".join(f"{k}={v}" for k, v in s["used"].items())})
    assert max_err < 1e-8, f"step energies differ from the research repo by {max_err:.2e} eV"

    out = HERE / "output"
    with open(out / f"{SURF}_steps.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    lines = [f"# Step free energies on {SURF} (eV), pH {SETTINGS['pH']}", "",
             "Written by `surface_steps.py`. dG = dG_raw + correction. Configurations and calc dirs: "
             f"`{SURF}_steps.csv`.", "",
             "| step | equation | U | dG_raw | correction | dG | correction terms |",
             "|---|---|---|---:|---:|---:|---|"]
    for r in rows:
        lines.append(f"| {r['step']} | {r['equation']} | {bias_label(r['bias'])} | {r['dG_raw_eV']:+.3f} | "
                     f"{r['correction_eV']:+.3f} | {r['dG_eV']:+.3f} | {r['correction_terms']} |")
    (out / f"{SURF}_steps.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[4:]))
    print(f"\nmatches the research step functions: max |diff| = {max_err:.1e} eV")
    print(f"wrote output/{SURF}_steps.csv, output/{SURF}_steps.md")


if __name__ == "__main__":
    main()
