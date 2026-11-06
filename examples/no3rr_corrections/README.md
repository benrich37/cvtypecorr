# NO3RR corrections: how they are fitted and applied (NO3_ads_gen)

This directory reproduces, using cvtypecorr only, the empirical corrections and the step free energies of the
NO3_ads_gen project's current nitrate-reduction pathway analysis (single-atom MN4-i sites; research repo,
`NO3_ads_gen/anl/mnc/no3rr_pathway_map.py`). It is meant for checking a different correction scheme against this one:

1. `fit_corrections.py` fits the species corrections (the chain in `corrections.get_collection`).
2. `surface_steps.py` computes the step free energies on one surface (FeN4-i), raw and corrected.
3. `no3_adsorption.py` gives NO3- adsorption free energies on all 29 surfaces: the comparison set.

Each script checks its results against the values the research repo's own code produced from the same data
(`data/research_reference.json`) and stops if they differ by more than 1e-8 eV. Currently all agree to 2e-12 eV.

```bash
cd examples/no3rr_corrections
python fit_corrections.py   # -> output/corrections.json
python surface_steps.py     # -> output/FeN4-i_steps.{csv,md}
python no3_adsorption.py    # -> output/no3_adsorption_energies.{csv,md}, output/no3_adsorption_all_configs.csv
```

Requirements: cvtypecorr and its dependencies (pymatgen, numpy, scipy). No research repo, no drive.

## Files

| file | contents |
|---|---|
| `energies.py` | free energy of a cache entry, molecule and adsorbed-state lookups, the step-energy rule (~100 lines) |
| `fit_corrections.py`, `surface_steps.py`, `no3_adsorption.py` | the three steps above |
| `data/settings.json` | T, kB, v0, pH, energy tags, electrode mu per bias, concentrations, provenance |
| `data/molecules.json` | the 11 free molecules used (cache entries + calc dirs) |
| `data/structures.json` | clean slabs (all surfaces), all 769 cached *NO3 configurations, and all configurations of *NO3H, *O, *OH, *NO2, *H on FeN4-i, at 0 and -1 V |
| `data/research_reference.json` | the research code's corrections, FeN4-i step energies and NO3 adsorption energies |
| `extract/extract_from_research.py` | how `data/` was made (needs the research repo and the drive) |
| `output/` | results of the three scripts |

## Method

### Energies (`energies.py`)

Every species, slab or molecule, is evaluated by its grand free energy at the electrode's electron chemical potential:

    G = E - TS - N*mu,   mu = -(U_SHE + 4.66 eV)

summing the cache terms `elec`, `trans` and `rot`. `elec` is the JDFTx energy and electronic smearing entropy of the
fixed-potential (slabs) or fixed-charge (molecules) calculation; `trans`/`rot` are the translational and rotational terms of free molecules,
solvated (or gas phase for `*_gas` molecules). **No vibrational or zero-point terms are used, for any species**, although
the cache has them for molecules. Because the electrons are counted in `-N*mu`, they never appear in a step: a step
written with `H3O+ -> H2O` is a proton-coupled electron transfer at potential U.

- Free molecules at a concentration c get `+kT ln c`: in the steps H3O+ at 10^-pH M (pH 6.3) and H2O at 55.34 M.
  NO3- and NO2- are at 1 M, H2 at 1 bar. T = 300 K.
- An adsorbed state `*X` is its lowest-G configuration among all cached ones, whatever its label.
- The fit uses mu = 0 (`No_bias`): electronless energies, compared with experiment at standard state.

### Fitting the corrections (`fit_corrections.py`)

A correction c(X) is added to the DFT free energy of species X. cvtypecorr's rigid fit fixes one species per
reaction so that `dE_DFT + sum(corrections) = dE_experiment`, using species already corrected. H2, H2O and HNO3
are grounded at 0. The DFT column is the raw DFT value of each reaction (before any correction):

| species | reaction | experiment | DFT | c (eV) |
|---|---|---:|---:|---:|
| H3O+ | 2 H3O+ + 2 e- -> H2 + 2 H2O | 0 V_SHE | -0.291 V | +0.2915 |
| NO3- | HNO3 + H2O -> NO3- + H3O+ | pKa -1.30 | -2.56 | -0.2171 |
| NO | NO3- + 4 H3O+ + 3 e- -> NO + 6 H2O | 0.960 V | 0.453 V | -0.5737 |
| NO2H | HNO2 + H3O+ + e- -> NO + 2 H2O | 0.996 V | 0.338 V | -0.2069 |
| NO2- | HNO2 + H2O -> NO2- + H3O+ | pKa 3.29 | 2.83 | -0.4710 |
| OH- | 2 H2O -> OH- + H3O+ | pKa 14.0 | 17.99 | -0.5274 |
| NH3 | NO2- + 7 H3O+ + 6 e- -> NH3 + 9 H2O | 0.792 V | 0.543 V | +0.0757 |
| NH4+ | NH4+ + H2O -> NH3 + H3O+ | pKa 9.24 | 19.51 | +0.9747 |

Placeholders, set to exactly 0 because no fit exists yet: O, H, N, NH, NH2, NOH, NHO, NH2O, NH3O, NHOH, NH2OH.

The H3O+ link aligns the proton with the SHE: after correction, the proton reference equals the computational
hydrogen electrode (1/2 H2 <-> H+ + e- at 0 V_SHE, pH 0). For *H on FeN4-i the two agree to 0.075 meV, which is
kT ln(55.5/55.34) from the water concentration (see the notes below).

### Applying them to a step (`surface_steps.py`)

Each adsorbate is labelled by the free species whose DFT error it is taken to share:
*NO3 -> NO3-, *NO2 -> NO2-, *OH -> OH-, *NO3H -> NO3H, *O -> O, *H -> H. The step correction is the sum of c(label)
over the products minus over the reactants, exactly as for a reaction between free species:

    correction = sum_products c(label) - sum_reactants c(label);   dG = dG_raw + correction

The step map of the research code passes these terms as `rct_tuples`/`pdt_tuples`, with the proton terms built into
each step function. Writing every species with its label gives the same numbers with no special cases. Three FeN4-i
steps at 0 V_SHE (all 20 are in `output/FeN4-i_steps.md`):

| step | dG_raw | correction | dG |
|---|---:|---:|---:|
| * + NO3- -> *NO3 | +0.193 | -c(NO3-) + c(NO3-) = 0 | +0.193 |
| * + NO3- -> *O + NO2- | +0.054 | -c(NO3-) + c(O) + c(NO2-) = -0.254 | -0.200 |
| *O + H3O+ -> *OH + H2O | -0.240 | -c(O) - c(H3O+) + c(OH-) + c(H2O) = -0.819 | -1.059 |

Adsorption and desorption steps therefore carry no correction. N-O cleavage and the protonation steps carry the
corrections of the free species they release and of the label change of the adsorbate. The last row shows what the
CHE alone would give (-0.291) and what the adsorbate labels add (*O -> *OH: c(OH-) - c(O) = -0.527 eV).

### NO3 adsorption for comparison (`no3_adsorption.py`)

`output/no3_adsorption_energies.md` (readable on its own) and `.csv` give, for each surface at 0 and -1 V_SHE:
dG_ads = G(*NO3) - G(*) - G(NO3-), the electronic-only dE_ads_elec, the electrons drawn from the electrode, the
configuration used with its label and convergence flag, the lowest intact and converged alternative, the three G
components, and the calc dirs. `no3_adsorption_all_configs.csv` has every configuration.

## Tracing a number back to the calculation

Every record in `data/*.json` carries:

- `calc_dir`, relative to the project root `PBED30_GBRV_20Ha`. The out file is `<calc_dir>/ion_opt/jdftx_run/out`.
  The root is `/Volumes/Ext4TB/scratch_backup/perl/NO3_ads_gen/PBED30_GBRV_20Ha` on the user's drive and
  `anvil:/anvil/projects/x-chm260113/Data/NO3_ads_gen/PBED30_GBRV_20Ha` on the group's Anvil copy.
- `out_file` size and mtime, to compare with the Anvil manifest (`PBED30_GBRV_20Ha.manifest.tsv.gz`, next to the
  database).
- `cache_key` and `cache_entry`: the key path in the research repo's energy cache `NO3_ads_gen/anl/mnc/data.json`
  and that entry, copied verbatim except that NaN became `null`.
- `label` (from `~/ht_data_anl/ht_labels.json`) and `ion_opt_finished` (`ion_opt/finished.txt`, or
  `sophie_calc_finished.txt` for calcs marked `sophie_mark.txt`).

`data/settings.json` -> `provenance` records the research commit (`6e77a2a`, code files unmodified), the sha1 of
`data.json` and of the label file, and the extraction time. For example, FeN4-i *NO3 at 0 V uses
`ads/FeN4-i/NO3/0V/in_bi_metal_atop_cc_bridge`, i.e.
`PBED30_GBRV_20Ha/ads/FeN4-i/NO3/0V/in_bi_metal_atop_cc_bridge/ion_opt/jdftx_run/out`.

To refresh `data/` after the cache changes, run `python extract/extract_from_research.py` in the research env with
the drive mounted, then rerun the three scripts.

## Notes for a comparison with another scheme

Points where two reasonable schemes are likely to differ, largest first:

1. **No vibrations or zero-point energy anywhere.** Including `vib` for the molecules alone would change the fitted
   corrections by up to 1.1 eV (H3O+ +0.093, NO3- -0.033, NO -0.741, NO2H -0.347, NO2- -0.426, OH- -0.335,
   NH3 -0.768, NH4+ -0.139). The adsorbates have no vibrational data in the cache.
2. **Translational and rotational free energy of the free molecule, none for the adsorbate.** Adsorbing NO3- costs
   its solvated trans+rot free energy, +0.508 eV at 300 K. This is why non-bonded (physisorbed) NO3 comes out at
   about +0.5 eV on the late and coinage metals, with dE_ads_elec close to 0.
3. **Adsorbates inherit the correction of a free species** (labels above). For NO3 adsorption the correction cancels.
   A scheme that corrects free NO3- but not *NO3 would get dG_ads - c(NO3-) = dG_ads + 0.217 eV. *O and *H have
   placeholder corrections of 0.
4. **Configuration choice.** The pathway analysis takes the lowest cached configuration whatever its label. On Ni,
   Cu, Pd, Ag, Ir, Pt, Au, Hg (and Co, Rh, Os at -1 V) that is a physisorbed or desorbed NO3, and a few minima come
   from unfinished optimizations. The `intact+fin` column gives the alternative. HgN4-i is in the table, but the
   pathway analysis leaves it out.
5. **D3.** The production JDFTx builds have the D3 C6 ordering error. Measured shifts on adsorption energies are about
   -7 meV per adsorbed O and +1 meV per N (research `Claude_Explore/05`), so roughly -20 meV for *NO3 (additive
   estimate), i.e. more exergonic than a corrected build would give. A per-calc assessment is in progress (research
   `NO3_ads_gen/anl/d3_c6_fix_impact/`, 2026-11-06).
6. **pH.** The fit is at standard state, so the corrections do not depend on pH: only H2O gets an activity there
   (55.5 M); H3O+ and the anions are at 1 M. pH enters the steps only through kT ln 10^-pH for H3O+.
7. **Temperature.** The cache is at 300 K, the experimental pKa's are at 298.15 K. In the acid links the DFT dG is
   converted to a pKa at 300 K and back to a dG at 298.15 K (as the research code does), which scales it by
   298.15/300. Making both 298.15 K moves the corrections by at most 7 meV (OH- -0.0066, NH4+ +0.0071, NO3-, NO,
   NO2H +0.0009).
8. **Water concentration.** 55.5 M in the fit, 55.34 M in the steps: a 0.075 meV difference.
9. **Charge drawn on adsorption.** On Ti, Zr, Hf and Nb at -1 V the lowest *NO3 configuration draws more than one extra
   electron from the electrode (dN_e +0.8 to +1.4). That is unusual (elsewhere -1.20 to +0.11). Check these structures
   before comparing them.
