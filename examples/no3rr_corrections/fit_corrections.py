"""Step 1: fit the species corrections used by the NO3RR pathway analysis.

Same chain, references and inputs as research.projects.no3_ads_gen.corrections.get_collection (the version used by
no3rr_pathway_map.py), written out with cvtypecorr only. Each link fixes the correction of one species from one
experimental number, using species whose corrections are already known:

    grounded at 0:  H2, H2O, NO3H
    H3O+  2 H3O+ + 2 e- -> H2 + 2 H2O        E0 = 0     V   (SHE)
    NO3-  pKa(HNO3) = -1.3
    NO    NO3- + 4 H3O+ + 3 e- -> NO + 6 H2O  E0 = 0.96  V
    NO2H  NO2H + H3O+ + e- -> NO + 2 H2O      E0 = 0.996 V
    NO2-  pKa(HNO2) = 3.29
    OH-   pKa(H2O)  = 14.0
    NH3   NO2- + 7 H3O+ + 6 e- -> NH3 + 9 H2O E0 = 0.792 V
    NH4+  pKa(NH4+) = 9.24
    placeholders at exactly 0 (no fit yet): O, H, N, NH, NH2, NOH, NHO, NH2O, NH3O, NHOH, NH2OH

DFT reaction energies use mu = 0 ("No_bias", electronless) and give only H2O an activity (55.5 M); see README.
A correction c(X) is added to the DFT free energy of X. The fit makes  dE_DFT + sum c = dE_exp  for each reaction.

Run: python fit_corrections.py   -> output/corrections.json, and checks against data/research_reference.json
"""
import json

from cvtypecorr.data.pka import data_dict, delta_g_to_pka
from cvtypecorr.data.rxn import AcidDeprotReaction, RedoxReaction
from cvtypecorr.fit.rigid import Collection

from energies import HERE, MOLECULES, REFERENCE, SETTINGS, V0, molecule_G, side_str

FIT_BIAS = "No_bias"
FIT_CONC = SETTINGS["conc_M"]["fit"]
GROUNDED = ["H2", "H2O", "NO3H"]
PLACEHOLDER_ZERO = ["O", "H", "N", "NH", "NH2", "NOH", "NHO", "NH2O", "NH3O", "NHOH", "NH2OH"]

# (fitted species, reaction type, reactants, products, experimental value: E0 in V vs SHE, or pKa)
CHAIN = [
    ("H3O+", "redox", [("H3O+", 2), ("e-", 2)], [("H2", 1), ("H2O", 2)], 0.0),
    ("NO3-", "acid", [("NO3H", 1), ("H2O", 1)], [("NO3-", 1), ("H3O+", 1)], data_dict["NO3H"]),
    ("NO", "redox", [("NO3-", 1), ("e-", 3), ("H3O+", 4)], [("NO", 1), ("H2O", 6)], 0.96),
    ("NO2H", "redox", [("NO2H", 1), ("e-", 1), ("H3O+", 1)], [("NO", 1), ("H2O", 2)], 0.996),
    ("NO2-", "acid", [("NO2H", 1), ("H2O", 1)], [("NO2-", 1), ("H3O+", 1)], data_dict["NO2H"]),
    ("OH-", "acid", [("H2O", 1), ("H2O", 1)], [("OH-", 1), ("H3O+", 1)], data_dict["H2O"]),
    ("NH3", "redox", [("NO2-", 1), ("e-", 6), ("H3O+", 7)], [("NH3", 1), ("H2O", 9)], 0.792),
    ("NH4+", "acid", [("NH4+", 1), ("H2O", 1)], [("NH3", 1), ("H3O+", 1)], data_dict["NH4+"]),
]


def dft_reaction(reactants, products) -> tuple[float, float]:
    """(dE, dN): electronless DFT free energy change and change in electron count, electrons themselves skipped."""
    dE = dN = 0.0
    for sign, side in ((-1, reactants), (+1, products)):
        for species, n in side:
            if species != "e-":
                dE += sign * n * molecule_G(species, FIT_BIAS, FIT_CONC.get(species))
                dN += sign * n * MOLECULES[species]["cache_entry"]["N"]
    return dE, dN


def fit_corrections(verbose: bool = False) -> tuple[Collection, list[dict]]:
    corrs = Collection(corrections={s: 0.0 for s in GROUNDED}, v0=V0)
    rows = []
    for species, kind, reactants, products, exp in CHAIN:
        dE, dN = dft_reaction(reactants, products)
        name = f"{side_str(reactants)} -> {side_str(products)}"
        if kind == "redox":
            dft = -dE / dN - V0  # equilibrium potential vs SHE that DFT predicts (dN = electrons consumed)
            corrs.apply_redox_correction(RedoxReaction(name, reactants, products, exp), dft_eq_pot=dft,
                                         overwrite=species)
        else:
            # Like the research code, the DFT dG is turned into a pKa at the cache temperature (300 K) and cvtypecorr
            # turns it back into a dG at the reaction's 298.15 K: the DFT dG enters scaled by 298.15/300.
            acid, base = reactants[0][0], products[0][0]
            dft = delta_g_to_pka(dE, n_protons=1, T=SETTINGS["T_K"])
            rxn = AcidDeprotReaction(name, acid, base, exp, ref_acid_base_pair_and_pKa=("H3O+", "H2O", 0.0),
                                     T=298.15)
            corrs.apply_acid_correction(rxn, dft, overwrite=species)
        rows.append({"species": species, "reaction": name, "kind": kind,
                     "exp": exp, "dft": dft, "unit": "V_SHE" if kind == "redox" else "pKa",
                     "correction": corrs.corrections[species]})
    for species in PLACEHOLDER_ZERO:
        corrs.corrections.setdefault(species, 0.0)
    if verbose:
        print(f"{'species':6s} {'experiment':>11s} {'DFT':>9s}  {'correction (eV)':>15s}  reaction")
        for r in rows:
            print(f"{r['species']:6s} {r['exp']:>7.3f} {r['unit']:<5s} {r['dft']:>7.3f}  {r['correction']:>+15.4f}  "
                  f"{r['reaction']}")
        print(f"grounded at 0: {', '.join(GROUNDED)};  placeholders at 0: {', '.join(PLACEHOLDER_ZERO)}")
    return corrs, rows


def check_against_research(corrs: Collection) -> float:
    ref = REFERENCE["corrections"]
    assert set(ref) == set(corrs.corrections), set(ref) ^ set(corrs.corrections)
    return max(abs(corrs.corrections[k] - ref[k]) for k in ref)


if __name__ == "__main__":
    corrs, rows = fit_corrections(verbose=True)
    err = check_against_research(corrs)
    assert err < 1e-8, f"corrections differ from the research repo by {err:.2e} eV"
    print(f"matches research corrections.get_collection: max |diff| = {err:.1e} eV")
    out = HERE / "output" / "corrections.json"
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps({
        "_about": "Species corrections (eV), added to the DFT free energy of each species. Written by "
                  "fit_corrections.py; 'fit' lists the reaction, the experimental value and the DFT value of each link.",
        "corrections": corrs.corrections, "v0_SHE_eV": corrs.v0, "fit": rows,
        "grounded_at_zero": GROUNDED, "placeholder_zero": PLACEHOLDER_ZERO}, indent=1) + "\n")
    print(f"wrote {out.relative_to(HERE)}")
