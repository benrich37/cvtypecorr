"""Step 2: free energy of each elementary step on one surface, before and after the corrections.

The steps are transitions among single-adsorbate states (*, *NO3, *NO3H, *O, *OH, *H, *NO2 through *NH3) in the step map of the NO3RR pathway
analysis (research: NO3_ads_gen/anl/mnc/no3rr_pathway_map.py), on FeN4-i at 0 and -1 V vs SHE, pH 6.3.

How a step energy is built (see energies.step_energy):
    dG_raw     = sum G(products) - sum G(reactants)
                 G(*X)  = lowest grand free energy of the cached *X configurations at the bias ("*" = clean slab)
                 G(mol) = grand free energy of the free molecule at the bias, + kT ln c for H3O+ (10^-pH M) and
                          H2O (55.34 M)
    correction = sum c(label) over products - over reactants, with each adsorbate labelled by the free species it
                 is taken to share its DFT error with (energies.adsorbate_label: *NO3 -> NO3-, *OH -> OH-, ...)
    dG         = dG_raw + correction
A proton is H3O+ -> H2O; its electron is implicit (all G are grand free energies at the electrode potential), so a
"H3O+" step is a proton-coupled electron transfer.

Run: python surface_steps.py   -> output/FeN4-i_steps.csv, output/FeN4-i_steps.md; checks against the research code
"""
import csv
import math

from energies import HERE, REFERENCE, SETTINGS, bias_label, side_str, step_energy
from fit_corrections import fit_corrections

SURF = REFERENCE["steps"]["surf"]
BIASES = ["0V", "n1p0V"]
CONC = SETTINGS["conc_M"]["steps"]

# (step map key, reactants, products). Every correction label of the step map appears in at least one step, except
# NH2O, NH2OH and NH3O: no adsorbed NH2O, NH2OH or NH3O has been calculated on any surface yet.
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
    ("NO2 -> NO2H", [("*NO2", 1), ("H3O+", 1)], [("*NO2H", 1), ("H2O", 1)]),
    ("NO2H -> NO", [("*NO2H", 1), ("H3O+", 1)], [("*NO", 1), ("H2O", 2)]),
    ("NO -> NHO", [("*NO", 1), ("H3O+", 1)], [("*NHO", 1), ("H2O", 1)]),
    ("NO -> NOH", [("*NO", 1), ("H3O+", 1)], [("*NOH", 1), ("H2O", 1)]),
    ("NHO -> NHOH", [("*NHO", 1), ("H3O+", 1)], [("*NHOH", 1), ("H2O", 1)]),
    ("NHOH -> NH", [("*NHOH", 1), ("H3O+", 1)], [("*NH", 1), ("H2O", 2)]),
    ("N -> NH", [("*N", 1), ("H3O+", 1)], [("*NH", 1), ("H2O", 1)]),
    ("NH -> NH2", [("*NH", 1), ("H3O+", 1)], [("*NH2", 1), ("H2O", 1)]),
    ("NH2 -> NH3", [("*NH2", 1), ("H3O+", 1)], [("*NH3", 1), ("H2O", 1)]),
]


def same(a: float | None, b: float | None) -> float:
    """|a - b|, with two missing values (NaN / null: a state without calcs) counting as equal."""
    missing = [v is None or math.isnan(v) for v in (a, b)]
    return 0.0 if all(missing) else (math.inf if any(missing) else abs(a - b))


def flags(used: dict) -> str:
    """Adsorbed states whose lowest configuration is not labelled intact, or whose ion_opt did not finish."""
    out = []
    for species, rec in used.items():
        if rec is None:
            out.append(f"{species}: no calc")
        elif species != "*" and (rec["label"].get("state") != "intact" or not rec["ion_opt_finished"]):
            state = rec["label"].get("state", "unlabelled")
            out.append(f"{species}: {state}{'' if rec['ion_opt_finished'] else ', unfinished'}")
    return "; ".join(out)


def main():
    corrections = fit_corrections()[0].corrections
    research = {(r["step"], r["bias"]): r for r in REFERENCE["steps"]["rows"]}
    rows, max_err = [], 0.0
    for key, reactants, products in STEPS:
        for bias in BIASES:
            s = step_energy(reactants, products, SURF, bias, corrections, CONC)
            ref = research[(key, bias)]
            max_err = max(max_err, same(s["dG_raw"], ref["dG_raw"]), same(s["dG"], ref["dG_corrected"]))
            rows.append({"step": key, "equation": f"{side_str(reactants)} -> {side_str(products)}", "bias": bias,
                         "U_vs_SHE_V": round(SETTINGS["bias_vs_SHE_V"][bias], 6),
                         "dG_raw_eV": s["dG_raw"], "correction_eV": s["correction"], "dG_eV": s["dG"],
                         "correction_terms": s["correction_terms"] or "0 (labels cancel)",
                         "flags": flags(s["used"]),
                         "calc_dirs": "; ".join(f"{k}={v['calc_dir'] if v else 'none'}" for k, v in s["used"].items())})
    assert max_err < 1e-8, f"step energies differ from the research repo by {max_err:.2e} eV"

    out = HERE / "output"
    with open(out / f"{SURF}_steps.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    lines = [f"# Step free energies on {SURF} (eV), pH {SETTINGS['pH']}", "",
             "Written by `surface_steps.py`. dG = dG_raw + correction. `flags` lists adsorbed states whose lowest "
             "configuration is not labelled intact or did not finish its optimization. Configurations and calc dirs: "
             f"`{SURF}_steps.csv`.", "",
             "| step | equation | U | dG_raw | correction | dG | correction terms | flags |",
             "|---|---|---|---:|---:|---:|---|---|"]
    fmt = lambda v: "n/a" if math.isnan(v) else f"{v:+.3f}"
    for r in rows:
        lines.append(f"| {r['step']} | {r['equation']} | {bias_label(r['bias'])} | {fmt(r['dG_raw_eV'])} | "
                     f"{r['correction_eV']:+.3f} | {fmt(r['dG_eV'])} | {r['correction_terms']} | {r['flags']} |")
    (out / f"{SURF}_steps.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[4:]))
    print(f"\nmatches the research step functions: max |diff| = {max_err:.1e} eV")
    print(f"wrote output/{SURF}_steps.csv, output/{SURF}_steps.md")


if __name__ == "__main__":
    main()
