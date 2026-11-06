"""Free energies of NO3_ads_gen calculations from the cache entries in data/ (standard library only).

This re-implements, in a few lines, what the research repo's `CachedNrgAnl` does with the settings of the current
NO3RR pathway analysis (research: NO3_ads_gen/anl/mnc/no3rr_pathway_map.py):

* Every species, slab or molecule, is evaluated by its grand free energy at the electrode's electron chemical
  potential mu:
        G = sum_t E[t] - sum_t TS[t] - N * mu        t in settings["energy_tags"] = elec, trans, rot
  so electrons never appear explicitly in a step: the electrons a step takes from the electrode are counted by -N*mu.
  "elec" is the DFT free energy (TS["elec"] is the smearing entropy); "trans"/"rot" are the solvated (or gas, for
  *_gas molecules) translational/rotational terms of free molecules. Vibrational terms ("vib") are in the cache but
  are NOT used, for any species.
* mu = settings["mu_eV"][bias]. "No_bias" means mu = 0 (plain, electronless free energies); the fit uses it.
* A free species at concentration c (M) gets + kT ln c. Species given no concentration are at 1 M (1 bar for H2).
* An adsorbed state is represented by its lowest-G configuration among all cached ones, with no filtering on labels,
  as in no3_pathway_sim_tools.get_lowest_struc_ref.
"""
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load(name: str) -> dict:
    return json.loads((HERE / "data" / f"{name}.json").read_text())


SETTINGS = load("settings")
MOLECULES = {r["species"]: r for r in load("molecules")["records"]}
STRUCTURES = load("structures")["records"]
REFERENCE = load("research_reference")

KT = SETTINGS["kB_eV_per_K"] * SETTINGS["T_K"]
TAGS = tuple(SETTINGS["energy_tags"])
MU = SETTINGS["mu_eV"]
V0 = SETTINGS["v0_SHE_eV"]

# Correction label of each adsorbate: the free species whose DFT error it is taken to share. Steps that only move a
# species between solution and surface (e.g. * + NO3- -> *NO3) therefore get no correction. Labels O and H are
# placeholders with correction exactly 0 (see README).
ADSORBATE_LABEL = {"NO3": "NO3-", "NO2": "NO2-", "OH": "OH-", "O": "O", "H": "H", "NO3H": "NO3H"}


def grand_free_energy(entry: dict, mu: float, tags=TAGS) -> float:
    """G = E - TS - N*mu of one verbatim cache entry, summing the given tags ("elec" must be present).

    NaN if a needed term is null (NaN in the cache; 4 entries here); such configurations are skipped, as in the
    research code.
    """
    if "elec" not in entry["E"]:
        raise KeyError("cache entry has no 'elec' energy")
    E = [entry["E"][t] for t in tags if t in entry["E"]]
    TS = [entry["TS"][t] for t in tags if t in entry["TS"]]
    if None in E + TS or entry["N"] is None:
        return math.nan
    return sum(E) - sum(TS) - entry["N"] * mu


def molecule_G(species: str, bias: str, conc: float | None = None, tags=TAGS) -> float:
    G = grand_free_energy(MOLECULES[species]["cache_entry"], MU[bias], tags)
    return G + (KT * math.log(conc) if conc else 0.0)


def configs(surf: str, ads: str, bias: str) -> list[dict]:
    """All records of adsorbate `ads` ("*" = clean slab) on `surf` at `bias`."""
    return [r for r in STRUCTURES if (r["surf"], r["ads"], r["bias"]) == (surf, ads, bias)]


def lowest(surf: str, ads: str, bias: str) -> tuple[float, dict]:
    """(G, record) of the lowest-G configuration of `ads` on `surf` at `bias`."""
    found = [(grand_free_energy(r["cache_entry"], MU[bias]), r) for r in configs(surf, ads, bias)]
    found = [(G, r) for G, r in found if math.isfinite(G)]
    if not found:
        raise LookupError(f"no {ads} configuration on {surf} at {bias}")
    return min(found, key=lambda pair: pair[0])


def step_energy(reactants, products, surf: str, bias: str, corrections: dict, conc: dict) -> dict:
    """Free energy of one step, raw and corrected.

    reactants/products: [(species, coefficient)], where "*X" is adsorbate X on `surf` ("*" the clean slab) and any
    other name is a free molecule (cvtypecorr name, e.g. "NO3-", "H3O+"). conc: {free species: concentration in M}.
    The correction is sum(products) - sum(reactants) of each species' correction label, i.e. what
    cvtypecorr's Collection._get_reaction_correction computes.
    """
    dG, corr, terms, used = 0.0, 0.0, [], {}
    for sign, side in ((-1, reactants), (+1, products)):
        for species, n in side:
            if species.startswith("*"):
                G, rec = lowest(surf, species[1:] or "*", bias)
                used[species] = rec["calc_dir"]
                label = ADSORBATE_LABEL.get(species[1:])
            else:
                G = molecule_G(species, bias, conc.get(species))
                label = species
            dG += sign * n * G
            if label is not None:
                corr += sign * n * corrections[label]
                terms.append(f"{'+' if sign > 0 else '-'}{n if n != 1 else ''}c({label})")
    return {"dG_raw": dG, "correction": corr, "dG": dG + corr, "correction_terms": " ".join(terms), "used": used}


def side_str(side) -> str:
    return " + ".join(s if n == 1 else f"{n} {s}" for s, n in side)


def bias_label(bias: str) -> str:
    return "No_bias" if bias == "No_bias" else f"{SETTINGS['bias_vs_SHE_V'][bias]:+.2f} V_SHE"
