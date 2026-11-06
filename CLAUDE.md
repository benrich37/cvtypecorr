# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Purpose

`cvtypecorr` derives per-species energy corrections for DFT calculations (CANDLE solvation by default) so that DFT-computed reaction energies reproduce experimental values — formation energies, standard reduction potentials, and pKa's.

## Commands

The package is installed editable in the `research` conda env (`/Users/richb/miniforge3/envs/research`), imported as `cvtypecorr` from `src/`.

```bash
pip install -e .                                    # install (editable)
pytest                                              # all tests
pytest tests/test_rigid.py::test_collection_acid    # single test
```

There is no linter/formatter config.

`/Users/richb/vs/research` (`src/research/projects/no3_ads_gen/corrections.py`, `get_collection()`) is a downstream consumer; changes to reaction energies or correction signs here change its fitted corrections.

## Architecture

Two layers: `data/` defines reactions with experimental reference energies; `fit/` fits species corrections against them.

### `cvtypecorr/data/rxn.py` — reactions
- Species are strings parsed by pymatgen `Composition`; charge is a `+`/`-` suffix (`"NO3-"`, `"H3O+"`, `"SO4-2"`). `"e-"` is the electron and is excluded from atom balance.
- Reactant/product lists accept `"H2O"` or `("H2O", 2)`; `resolve_specie_set` normalizes to `[(name, coeff), ...]`. The `Reaction` base class raises `ValueError` if it is not atom- and charge-balanced.
- Every reaction exposes `delta_E` (eV) as the reference value the fitter compares against:
  - `ThermoReaction` — `delta_E = delta_g` given directly.
  - `FormationReaction` — builds reactants from `reference_species` (O←H2O, H←H2, N←NH3, C←CH4) by iterating elements in `el_by_priority` (default: descending Z) until mass balance; reference species can introduce elements absent from the target formula (with possibly negative coefficients). Adding an element requires a `reference_species` entry.
  - `RedoxReaction` — electrons must appear explicitly (reactants for reduction). `delta_E` is the *electronless* reaction energy derived from the potential vs SHE and `v0` (absolute SHE potential, default `v0_SHE_CANDLE = 4.66`); setting `.v0` recomputes `delta_E`. `RedoxReaction.electrons` is the summed electron *charge* (reactants minus products), so it is **negative** for reductions; the conversion functions (`eq_pot_to_electronless_reduction_energy` and its inverse) expect that signed value, giving `delta_E = -|n|(V + v0)` for a reduction.
  - `AcidDeprotReaction` — `acid + ref_base -> base + ref_acid` (default reference pair H3O+/H2O with pKa 0.0); `delta_E` is the pKa difference converted to free energy at temperature `T`.
- `data/pka.py` holds `pka_to_delta_g`/`delta_g_to_pka` and a small table of reference pKa's (`data_dict`).

### `cvtypecorr/fit/rigid.py` — `Collection`
- `Collection.corrections` is a plain `dict[str, float]` of additive per-species energy corrections, plus `v0`.
- "Rigid" fitting: `apply_correction(rxn, dft_nrg)` solves for exactly **one** species missing from `corrections`; it raises if more than one is missing. Corrections must therefore be chained species-by-species, grounding some species (typically neutrals) to 0 first. `overwrite=<species>` drops that species' existing correction so it is refit.
- `dft_nrg` is the **raw, uncorrected** DFT reaction energy; already-known corrections are folded in via `get_reaction_correction` (products minus reactants, weighted by coefficient). The new correction is chosen so that `dft_nrg + get_reaction_correction(rxn) == rxn.delta_E`.
- `apply_redox_correction` accepts a DFT equilibrium potential or an electronless reaction energy (exactly one), syncs `rxn.v0` to the collection's `v0`, and pins `corrections["e-"] = 0`. `apply_acid_correction` takes a DFT pKa.

`examples/obtaining_rigid_corrections.ipynb` walks through a full correction chain (H3O+ → NO3- → NO → NO2H → NO2-) and is the best reference for intended usage and sign conventions.

`examples/no3rr_corrections/` (branch `no3rr-correction-example`, 2026-11-06) is the same chain on real data: it fits the NO3_ads_gen corrections, applies them to adsorption steps on FeN4-i and gives NO3- adsorption energies on 29 surfaces. It runs without the research repo, and its three scripts assert agreement with the research repo's own results (`data/research_reference.json`), so run them after changing reaction energies or correction signs here. Its README has the method, data provenance and the scheme-comparison notes; `extract/extract_from_research.py` regenerates `data/`.
