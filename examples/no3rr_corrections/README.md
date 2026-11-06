# NO3RR corrections: how they are fitted and applied (NO3_ads_gen)

This directory reproduces, using cvtypecorr only, the empirical corrections and the step free energies of the
NO3_ads_gen project's current nitrate-reduction pathway analysis (single-atom MN4-i sites; research repo,
`NO3_ads_gen/anl/mnc/no3rr_pathway_map.py`). It is meant for checking a different correction scheme against this one:

1. `fit_corrections.py` fits the species corrections (the chain in `corrections.get_collection`): every N
   intermediate of the pathway gets one, and closure checks test the result against couples not used in the fit.
2. `surface_steps.py` computes the free energies of 19 steps on one surface (FeN4-i), raw and corrected.
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
| `energies.py` | free energy of a cache entry, molecule and adsorbed-state lookups, the step-energy rule (~120 lines) |
| `fit_corrections.py`, `surface_steps.py`, `no3_adsorption.py` | the three steps above |
| `data/settings.json` | T, kB, v0, pH, energy tags, electrode mu per bias, concentrations, provenance |
| `data/molecules.json` | the 24 free molecules used (14 aqueous/reference calcs, 10 gas-phase twins `X(g)`), cache entries + calc dirs |
| `data/structures.json` | clean slabs (all surfaces), all 769 cached *NO3 configurations, and all configurations of the 15 adsorbates in the FeN4-i steps, at 0 and -1 V |
| `data/research_reference.json` | the research code's corrections, FeN4-i step energies and NO3 adsorption energies |
| `extract/extract_from_research.py` | how `data/` was made (needs the research repo and the drive) |
| `output/` | results of the three scripts |

## Method

### Energies (`energies.py`)

Every species, slab or molecule, is evaluated by its grand free energy at the electrode's electron chemical potential:

    G = E - TS - N*mu,   mu = -(U_SHE + 4.66 eV)

summing the cache terms `elec`, `trans` and `rot`. `elec` is the JDFTx energy and electronic smearing entropy of the
fixed-potential (slabs) or fixed-charge (molecules) calculation; `trans`/`rot` are the translational and rotational
terms of free molecules, solvated (or gas phase for `*_gas` molecules). **No vibrational or zero-point terms are
used, for any species**, although the cache has them for molecules. Because the electrons are counted in `-N*mu`,
they never appear in a step: a step written with `H3O+ -> H2O` is a proton-coupled electron transfer at potential U.

- Free molecules at a concentration c get `+kT ln c`: in the steps H3O+ at 10^-pH M (pH 6.3) and H2O at 55.34 M.
  NO3- and NO2- are at 1 M, H2 at 1 bar. T = 300 K.
- An adsorbed state `*X` is its lowest-G configuration among all cached ones, whatever its label.
- The fit uses mu = 0 (`No_bias`): electronless energies, compared with experiment at standard state.

### Fitting the corrections (`fit_corrections.py`)

A correction c(X) is added to the DFT free energy of species X. cvtypecorr's rigid fit fixes one species per
reaction so that `dE_DFT + sum(corrections) = dE_experiment`, using species already corrected. H2, H2O and HNO3
are grounded at 0. The DFT column is the raw DFT value of each reaction (before any correction).

The first eight links are the original chain. The rest are the NHxOy extension (research, 2026-11-06), fitted to the
benchmark report `NO3_ads_gen/gen/from_aziz/Experimental Benchmark Thermochemistry for NHxOyHz Species in
Nitrate-to-Ammonia Reduction_ Correcting DFT Gas-Phase Errors.pdf` ("the report", sections in brackets) with the
reference molecules of `gen/from_aziz/sub_missing_molecules.py`:

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
| NHO | NH3 + H2O -> HNO + 2 H2 (aq formation) [C] | dG 3.939 eV | 4.409 | -0.3936 |
| NH2O | NH3 + H2O -> H2NO + 1.5 H2 (aq formation) [C] | dG 3.416 eV | 3.372 | +0.1197 |
| NH2OH | H2NO + H3O+ + e- -> NH2OH + H2O [D] | 0.900 V | 0.497 V | +0.0086 |
| NH2 | NH3 -> NH2 + 0.5 H2 (gas) [A] | dH 2.400 eV | 2.634 | -0.1582 |
| NH | NH3 -> NH + H2 (gas) [A] | dH 4.191 eV | 4.680 | -0.4139 |
| N | NH3 -> N + 1.5 H2 (gas, quartet N) [A] | dH 5.369 eV | 6.226 | -0.7820 |
| NOH | HNO -> NOH (gas) [A] | dH 1.121 eV | 0.983 | -0.2554 |
| NHOH | H2NO -> HNOH (gas) [A] | dH 0.315 eV | 0.450 | -0.0152 |
| NH3O | NH2OH -> NH3O (gas) [A] | dH 1.077 eV | 0.893 | +0.1930 |

- HNO(aq) is NO(aq) (+102.0 kJ/mol, Armstrong 2015) plus E(NO,H+/HNO) = -0.15 V; H2NO(aq) is +66 kJ/mol
  (Armstrong 2015); NH3(aq) and H2O(l) are NBS values.
- The six gas-phase links compare ATcT 298 K enthalpies with the `*_gas` calcs (E terms plus pV). Each one moves
  the gas-phase error of a reaction from an already-corrected partner onto the new species, so the new
  correction is the partner's correction plus that error. For the six species with both aqueous and gas data, the
  two kinds of error agree to within 0.08 eV (research check of 2026-10-23). The ATcT values for NOH, HNOH and NH3O are theory-dominated, and an
  atom (N) is a poor GGA test.
- Zero corrections remain only for metal-bound O and H. They have no molecular analog, and H2O and H2 are the
  references.

The H3O+ link aligns the proton with the SHE: after correction, the proton reference equals the computational
hydrogen electrode (1/2 H2 <-> H+ + e- at 0 V_SHE, pH 0). For *H on FeN4-i the two agree to 0.075 meV, which is
kT ln(55.5/55.34) from the water concentration.

**Closure checks** (couples not used in the fit, recomputed with the corrected energies; `output/corrections.json`):

| couple | experiment (V) | corrected | raw DFT | difference |
|---|---:|---:|---:|---:|
| HNO,H+/H2NO (Armstrong 2015) [D] | 0.520 | 0.523 | 0.745 | +3 mV |
| NO3- + 2H+ + 2e- -> NO2- (from NBS dfG) [D] | 0.832 | 0.844 | 0.426 | +12 mV |
| NO2- + 8H+ + 6e- -> NH4+ (from NBS dfG) [D] | 0.901 | 0.884 | 0.737 | -17 mV |
| NO3- + 10H+ + 8e- -> NH4+ (from NBS dfG) [D] | 0.884 | 0.874 | 0.659 | -10 mV |

### Applying them to a step (`surface_steps.py`)

Each adsorbate is labelled by the free species whose DFT error it is taken to share. *NO3, *NO2 and *OH carry the
anion's correction (NO3-, NO2-, OH-). Every other adsorbate carries the correction of its own name (*NO -> NO,
*NHO -> NHO, ..., *NH3 -> NH3), and *O and *H carry 0. The step correction is the sum of c(label) over the products
minus over the reactants, exactly as for a reaction between free species:

    correction = sum_products c(label) - sum_reactants c(label);   dG = dG_raw + correction

The step map of the research code passes these terms as `rct_tuples`/`pdt_tuples`, with the proton terms built into
each step function. Writing every species with its label gives the same numbers with no special cases. Four FeN4-i
steps at 0 V_SHE (all 19 steps at both potentials are in `output/FeN4-i_steps.md`):

| step | dG_raw | correction | dG |
|---|---:|---:|---:|
| * + NO3- -> *NO3 | +0.193 | -c(NO3-) + c(NO3-) = 0 | +0.193 |
| * + NO3- -> *O + NO2- | +0.054 | -c(NO3-) + c(O) + c(NO2-) = -0.254 | -0.200 |
| *O + H3O+ -> *OH + H2O | -0.240 | -c(O) - c(H3O+) + c(OH-) + c(H2O) = -0.819 | -1.059 |
| *NO + H3O+ -> *NHO + H2O | +0.849 | -c(NO) - c(H3O+) + c(NHO) + c(H2O) = -0.111 | +0.738 |

Adsorption and desorption steps therefore carry no correction. N-O cleavage and the protonation steps carry the
corrections of the free species they release and of the change in the adsorbate's label. In the *O -> *OH row,
-0.291 eV is what the CHE alone would give (the H3O+ term); the adsorbate labels add c(OH-) - c(O) = -0.527 eV.

The steps cover every correction label of the step map except NH2O, NH2OH and NH3O: no adsorbed NH2O, NH2OH or NH3O
has been calculated on any surface. The `flags` column marks steps whose FeN4-i minimum is unfinished: *NO2H at 0 V,
*NHO and *NOH at -1 V, and *NH3. It also marks *NOH at 0 V, which has no calc.

### NO3 adsorption for comparison (`no3_adsorption.py`)

`output/no3_adsorption_energies.md` (readable on its own) and `.csv` give, for each surface at 0 and -1 V_SHE:
dG_ads = G(*NO3) - G(*) - G(NO3-), the electronic-only dE_ads_elec, the electrons drawn from the electrode, the
configuration used with its label and convergence flag, the lowest intact and converged alternative, the three G
components, and the calc dirs. `no3_adsorption_all_configs.csv` has every configuration. None of these numbers depend
on the NHxOy extension.

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

`data/settings.json` -> `provenance` records:
- the research commit (`6e77a2a`);
- the sha1 of each research file the numbers depend on. Their 2026-11-06 changes are not committed there yet: the
  extension in `corrections.py`, and the `data.json` re-cache of NH2O_gas, N_gas and N;
- the sha1 of the label file;
- the extraction time.

For example, FeN4-i *NO3 at 0 V uses `ads/FeN4-i/NO3/0V/in_bi_metal_atop_cc_bridge`, i.e.
`PBED30_GBRV_20Ha/ads/FeN4-i/NO3/0V/in_bi_metal_atop_cc_bridge/ion_opt/jdftx_run/out`.

To refresh `data/` after the cache changes, run `python extract/extract_from_research.py` in the research env with
the drive mounted, then rerun the three scripts.

## Notes for a comparison with another scheme

Points where two reasonable schemes are likely to differ, largest first:

1. **No vibrations or zero-point energy anywhere.** Fitting with `vib` included for the molecules would change the
   corrections by up to 1.1 eV (NH4+). The `full` values are H3O+ +0.093, NO3- -0.033, NO -0.741, NO2H -0.347,
   NO2- -0.426, OH- -0.335, NH3 -0.768, NH4+ -0.139, NHO -0.675, NH2O -0.348, NH2OH -0.684, NH2 -0.731,
   NH -0.814, N -1.118, NOH -0.547, NHOH -0.493 and NH3O -0.479. The adsorbates have no vibrational data in the
   cache, so the corrections in use absorb the zero-point energy that is left out.
2. **Translational and rotational free energy of the free molecule, none for the adsorbate.** Adsorbing NO3- costs
   its solvated trans+rot free energy, +0.508 eV at 300 K. This is why non-bonded (physisorbed) NO3 comes out at
   about +0.5 eV on the late and coinage metals, with dE_ads_elec close to 0.
3. **Adsorbates inherit the correction of a free species** (labels above). For NO3 adsorption the correction cancels.
   A scheme that corrects free NO3- but not *NO3 would get dG_ads - c(NO3-) = dG_ads + 0.217 eV.
4. **Configuration choice.** The pathway analysis takes the lowest cached configuration whatever its label. On Ni,
   Cu, Pd, Ag, Ir, Pt, Au, Hg (and Co, Rh, Os at -1 V) that is a physisorbed or desorbed NO3, and a few minima come
   from unfinished optimizations. The `intact+fin` column gives the alternative. HgN4-i is in the table, but the
   pathway analysis leaves it out.
5. **Reference data.** The chain is not built from one self-consistent data set:
   - The extension combines NBS dfG (NH3, H2O) with Armstrong values (HNO, H2NO, NH2OH). The report (section C)
     advises against mixing these two sources.
   - The original links use E(NO3-/NO) 0.96 V, E(HNO2/NO) 0.996 V, pKa(HNO2) 3.29 and pKa(NH4+) 9.24. The report
     gives 0.957 V, 0.983 V (sources conflict, 0.983-1.036), 3.15 and 9.25.
   - The NO2-/NH3 couple (0.792 V) comes from a separate source (doi:10.1002/elsa.202100142).
   - The closure checks put the resulting inconsistency at 10-17 mV per electron on the multi-electron NBS couples.
   - A rebuild on a single data set is a possible later revision: formation reactions against fixed anchors, as in
     the research repo's Patty project.
6. **D3.** All these calculations, and every production JDFTx build so far, have the D3 C6 indexing error (fixed
   in jdftx b758dd5c). The fix makes these NO3 adsorption energies less exergonic by +12 meV median (0 to +79 meV;
   largest on TcN4-i, +41 meV on FeN4-i at 0 V). That comes from the exact per-calc D3 energy corrections
   `dE_final` in research `NO3_ads_gen/anl/d3_c6_fix_impact/summary.csv`. 57 of the 58 configurations used need only
   that energy correction; CoN4-i at 0 V would need re-optimization. Energies from the same production builds are
   consistent with these; a fixed build or another code would differ by about this much.
7. **pH.** The fit is at standard state, so the corrections do not depend on pH: only H2O gets an activity there
   (55.5 M); H3O+ and the anions are at 1 M. pH enters the steps only through kT ln 10^-pH for H3O+.
8. **Temperature.** The cache is at 300 K, the experimental data at 298.15 K. In the acid links the DFT dG is
   converted to a pKa at 300 K and back to a dG at 298.15 K (as the research code does), which scales it by
   298.15/300. Making both 298.15 K moves the corrections by at most 7 meV (OH- -0.0066, NH4+ +0.0071; NO3-, NO
   and NO2H +0.0009).
9. **Water concentration.** 55.5 M in the fit, 55.34 M in the steps: a 0.075 meV difference.
10. **Charge drawn on adsorption.** On Ti, Zr, Hf and Nb at -1 V the lowest *NO3 configuration draws more than one
    extra electron from the electrode (dN_e +0.8 to +1.4). That is unusual (elsewhere -1.20 to +0.11). Check these
    structures before comparing them.
