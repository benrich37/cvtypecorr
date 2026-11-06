"""Step 1: fit the species corrections used by the NO3RR pathway analysis.

Same chain, references and inputs as research.projects.no3_ads_gen.corrections.get_collection (the version used by
no3rr_pathway_map.py), written out with cvtypecorr only. Each link fixes the correction of one species from one
experimental number, using species whose corrections are already known:

    grounded at 0:  H2, H2O, NO3H
    H3O+   2 H3O+ + 2 e- -> H2 + 2 H2O         E0 = 0     V   (SHE)
    NO3-   pKa(HNO3) = -1.3
    NO     NO3- + 4 H3O+ + 3 e- -> NO + 6 H2O   E0 = 0.96  V
    NO2H   NO2H + H3O+ + e- -> NO + 2 H2O       E0 = 0.996 V
    NO2-   pKa(HNO2) = 3.29
    OH-    pKa(H2O)  = 14.0
    NH3    NO2- + 7 H3O+ + 6 e- -> NH3 + 9 H2O  E0 = 0.792 V
    NH4+   pKa(NH4+) = 9.24
    NHxOy extension (research, 2026-11-06), every other N intermediate of the pathway:
    NHO, NH2O         aqueous formation from NH3 + H2O                dfG(aq)
    NH2OH  NH2O + H3O+ + e- -> NH2OH + H2O      E0 = 0.900 V
    NH2, NH, N        gas-phase NH3 -> X + n H2                       ATcT dfH(gas)
    NOH, NHOH, NH3O   gas-phase isomerization of NHO, NH2O, NH2OH     ATcT dfH(gas), theory-dominated
    zero (no molecular analog): metal-bound O and H

DFT reaction energies use mu = 0 ("No_bias", electronless) and give only H2O an activity (55.5 M); the gas-phase links
compare enthalpies (E terms of the *_gas calcs plus pV), see README. A correction c(X) is added to the DFT free energy
of X. The fit makes  dE_DFT + sum c = dE_exp  for each reaction.

Run: python fit_corrections.py   -> output/corrections.json, and checks against data/research_reference.json
"""
import json

from cvtypecorr.data.pka import data_dict, delta_g_to_pka
from cvtypecorr.data.rxn import AcidDeprotReaction, RedoxReaction, ThermoReaction
from cvtypecorr.fit.rigid import Collection

from energies import HERE, KT, MOLECULES, REFERENCE, SETTINGS, TAGS, V0, molecule_G, side_str

FIT_BIAS = "No_bias"
FIT_CONC = SETTINGS["conc_M"]["fit"]
GROUNDED = ["H2", "H2O", "NO3H"]
ZERO = ["O", "H"]

# Reference thermochemistry of the NHxOy links (kJ/mol, 298.15 K), as in research corrections.py. Source: the
# benchmark report research:NO3_ads_gen/gen/from_aziz/"Experimental Benchmark Thermochemistry for NHxOyHz Species in
# Nitrate-to-Ammonia Reduction_ Correcting DFT Gas-Phase Errors.pdf", section in brackets.
KJ_PER_EV = 96.485
DFG_AQ = {"H2": 0.0, "NH3": -26.50, "H2O": -237.129,  # NBS (H2O liquid) [C]
          "NHO": 102.0 + 0.15 * KJ_PER_EV,    # HNO(aq): NO(aq) +102.0 (Armstrong 2015) and E(NO,H+/HNO) = -0.15 V [C, D]
          "NH2O": 66.0}                       # H2NO(aq), Armstrong 2015, +/- 3 [C]
DFH_GAS = {"H2": 0.0, "NH3": -45.556, "NH2": 186.03, "NH": 358.80, "N": 472.445,  # ATcT v1.202 [A]
           "NHO": 106.96, "NOH": 215.12,      # NOH: v1.122e, theory-dominated
           "NH2O": 64.05, "NHOH": 94.43,      # trans-HNOH, theory-dominated
           "NH2OH": -43.44, "NH3O": 60.5}     # trans-NH2OH; NH3O 59.9-61.1 (v1.122) midpoint, theory-dominated


def table_value(table, reactants, products) -> float:
    """Reaction value in eV from a table of formation values in kJ/mol."""
    return (sum(n * table[s] for s, n in products) - sum(n * table[s] for s, n in reactants)) / KJ_PER_EV


# (fitted species, reaction type, reactants, products, experimental value: E0 in V vs SHE, pKa, or dG/dH in eV)
CHAIN = [
    ("H3O+", "redox", [("H3O+", 2), ("e-", 2)], [("H2", 1), ("H2O", 2)], 0.0),
    ("NO3-", "acid", [("NO3H", 1), ("H2O", 1)], [("NO3-", 1), ("H3O+", 1)], data_dict["NO3H"]),
    ("NO", "redox", [("NO3-", 1), ("e-", 3), ("H3O+", 4)], [("NO", 1), ("H2O", 6)], 0.96),
    ("NO2H", "redox", [("NO2H", 1), ("e-", 1), ("H3O+", 1)], [("NO", 1), ("H2O", 2)], 0.996),
    ("NO2-", "acid", [("NO2H", 1), ("H2O", 1)], [("NO2-", 1), ("H3O+", 1)], data_dict["NO2H"]),
    ("OH-", "acid", [("H2O", 1), ("H2O", 1)], [("OH-", 1), ("H3O+", 1)], data_dict["H2O"]),
    ("NH3", "redox", [("NO2-", 1), ("e-", 6), ("H3O+", 7)], [("NH3", 1), ("H2O", 9)], 0.792),
    ("NH4+", "acid", [("NH4+", 1), ("H2O", 1)], [("NH3", 1), ("H3O+", 1)], data_dict["NH4+"]),
    *[(x, "formation", r, p, table_value(DFG_AQ, r, p)) for x, r, p in [
        ("NHO", [("NH3", 1), ("H2O", 1)], [("NHO", 1), ("H2", 2)]),
        ("NH2O", [("NH3", 1), ("H2O", 1)], [("NH2O", 1), ("H2", 1.5)])]],
    ("NH2OH", "redox", [("NH2O", 1), ("e-", 1), ("H3O+", 1)], [("NH2OH", 1), ("H2O", 1)], 0.900),
    *[(x, "gas", r, p, table_value(DFH_GAS, r, p)) for x, r, p in [
        ("NH2", [("NH3", 1)], [("NH2", 1), ("H2", 0.5)]),
        ("NH", [("NH3", 1)], [("NH", 1), ("H2", 1)]),
        ("N", [("NH3", 1)], [("N", 1), ("H2", 1.5)]),
        ("NOH", [("NHO", 1)], [("NOH", 1)]),
        ("NHOH", [("NH2O", 1)], [("NHOH", 1)]),
        ("NH3O", [("NH2OH", 1)], [("NH3O", 1)])]],
]

# Couples not used in the fit, recomputed with the corrected energies (E0 in V vs SHE, report section)
CLOSURE = [
    ("HNO,H+/H2NO", [("NHO", 1), ("H3O+", 1), ("e-", 1)], [("NH2O", 1), ("H2O", 1)], 0.52, "Armstrong 2015 [D]"),
    ("NO3-/NO2-", [("NO3-", 1), ("H3O+", 2), ("e-", 2)], [("NO2-", 1), ("H2O", 3)], 0.832, "from NBS dfG [D]"),
    ("NO2-/NH4+", [("NO2-", 1), ("H3O+", 8), ("e-", 6)], [("NH4+", 1), ("H2O", 10)], 0.901, "from NBS dfG [D]"),
    ("NO3-/NH4+", [("NO3-", 1), ("H3O+", 10), ("e-", 8)], [("NH4+", 1), ("H2O", 13)], 0.884, "from NBS dfG [D]"),
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


def gas_enthalpy_change(reactants, products) -> float:
    """dH of a gas-phase reaction from the *_gas calcs ("X(g)" records; H2 is gas already): E terms + kT per mol gas."""
    def h(species):
        entry = MOLECULES["H2" if species == "H2" else f"{species}(g)"]["cache_entry"]
        return sum(entry["E"][t] for t in TAGS if t in entry["E"])
    dH = sum(n * h(s) for s, n in products) - sum(n * h(s) for s, n in reactants)
    return dH + KT * (sum(n for _, n in products) - sum(n for _, n in reactants))


def fit_corrections(verbose: bool = False) -> tuple[Collection, list[dict]]:
    corrs = Collection(corrections={s: 0.0 for s in GROUNDED}, v0=V0)
    rows = []
    for species, kind, reactants, products, exp in CHAIN:
        name = f"{side_str(reactants)} -> {side_str(products)}"
        if kind == "redox":
            dE, dN = dft_reaction(reactants, products)
            dft = -dE / dN - V0  # equilibrium potential vs SHE that DFT predicts (dN = electrons consumed)
            corrs.apply_redox_correction(RedoxReaction(name, reactants, products, exp), dft_eq_pot=dft,
                                         overwrite=species)
        elif kind == "acid":
            # Like the research code, the DFT dG is turned into a pKa at the cache temperature (300 K) and cvtypecorr
            # turns it back into a dG at the reaction's 298.15 K: the DFT dG enters scaled by 298.15/300.
            acid, base = reactants[0][0], products[0][0]
            dft = delta_g_to_pka(dft_reaction(reactants, products)[0], n_protons=1, T=SETTINGS["T_K"])
            rxn = AcidDeprotReaction(name, acid, base, exp, ref_acid_base_pair_and_pKa=("H3O+", "H2O", 0.0),
                                     T=298.15)
            corrs.apply_acid_correction(rxn, dft, overwrite=species)
        else:
            dft = dft_reaction(reactants, products)[0] if kind == "formation" else gas_enthalpy_change(reactants, products)
            corrs.apply_correction(ThermoReaction(name, reactants, products, exp), dft, overwrite=species)
        unit = {"redox": "V_SHE", "acid": "pKa", "formation": "eV dG", "gas": "eV dH(g)"}[kind]
        rows.append({"species": species, "reaction": name, "kind": kind, "exp": exp, "dft": dft, "unit": unit,
                     "correction": corrs.corrections[species]})
    for species in ZERO:
        corrs.corrections.setdefault(species, 0.0)
    if verbose:
        print(f"{'species':7s} {'experiment':>14s} {'DFT':>8s}  {'c (eV)':>8s}  reaction")
        for r in rows:
            print(f"{r['species']:7s} {r['exp']:>7.3f} {r['unit']:<8s} {r['dft']:>6.3f}  {r['correction']:>+8.4f}  "
                  f"{r['reaction']}")
        print(f"grounded at 0: {', '.join(GROUNDED)};  zero (no molecular analog): {', '.join(ZERO)}")
    return corrs, rows


def closure_checks(corrs: Collection) -> list[dict]:
    """Corrected and raw DFT equilibrium potentials of the couples in CLOSURE."""
    rows = []
    for couple, reactants, products, exp, source in CLOSURE:
        dE, dN = dft_reaction(reactants, products)
        c = (sum(n * corrs.corrections[s] for s, n in products)
             - sum(n * corrs.corrections[s] for s, n in reactants if s != "e-"))
        rows.append({"couple": couple, "reaction": f"{side_str(reactants)} -> {side_str(products)}", "exp_V": exp,
                     "corrected_V": -(dE + c) / dN - V0, "raw_V": -dE / dN - V0, "source": source})
    return rows


def check_against_research(corrs: Collection) -> float:
    ref = REFERENCE["corrections"]
    assert set(ref) == set(corrs.corrections), set(ref) ^ set(corrs.corrections)
    return max(abs(corrs.corrections[k] - ref[k]) for k in ref)


if __name__ == "__main__":
    corrs, rows = fit_corrections(verbose=True)
    err = check_against_research(corrs)
    assert err < 1e-8, f"corrections differ from the research repo by {err:.2e} eV"
    print(f"matches research corrections.get_collection: max |diff| = {err:.1e} eV\n")
    closure = closure_checks(corrs)
    print(f"{'closure check':13s} {'exp (V)':>8s} {'corrected':>10s} {'raw':>7s}  {'diff (mV)':>9s}  source")
    for r in closure:
        print(f"{r['couple']:13s} {r['exp_V']:>8.3f} {r['corrected_V']:>10.3f} {r['raw_V']:>7.3f}  "
              f"{1000 * (r['corrected_V'] - r['exp_V']):>+9.0f}  {r['source']}")
    out = HERE / "output" / "corrections.json"
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps({
        "_about": "Species corrections (eV), added to the DFT free energy of each species. Written by "
                  "fit_corrections.py; 'fit' lists the reaction, the experimental value and the DFT value of each link, "
                  "'closure_checks' couples not used in the fit.",
        "corrections": corrs.corrections, "v0_SHE_eV": corrs.v0, "fit": rows, "closure_checks": closure,
        "grounded_at_zero": GROUNDED, "zero_no_molecular_analog": ZERO}, indent=1) + "\n")
    print(f"wrote {out.relative_to(HERE)}")
