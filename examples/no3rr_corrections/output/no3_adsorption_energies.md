# NO3- adsorption free energies on MN4-i single-atom sites

Reaction `* + NO3-(aq) -> *NO3` at fixed electrode potential, from the NO3_ads_gen calculations (PBE-D3(0), GBRV, 20 Ha, CANDLE LinearPCM, 0.1 M NaF, grand canonical at fixed U). Energies in eV. Written by `no3_adsorption.py`; every number is in `no3_adsorption_energies.csv` with its components and calc dirs (`calc_dir`, relative to `PBED30_GBRV_20Ha`, on Anvil at `anvil:/anvil/projects/x-chm260113/Data/NO3_ads_gen/PBED30_GBRV_20Ha`), and every configuration in `no3_adsorption_all_configs.csv`.

**Definition.** dG_ads = G(*NO3) - G(*) - G(NO3-), each G a grand free energy at the electrode potential, G = E - TS - N*mu with mu = -(U_SHE + 4.66) eV, so the electrons drawn from the electrode (dN_e) are counted in G and no separate electron term appears.

- G(*), G(*NO3): DFT free energy of the slab (E - TS of the electronic smearing). No vibrations or zero-point energy. G(*NO3) is the lowest of all cached *NO3 configurations.
- G(NO3-): DFT free energy of solvated NO3- plus its solvated translational and rotational free energy (-0.508 eV at 300 K), at 1 M. No vibrations or zero-point energy.
- Correction: none in effect. NO3- carries an empirical correction (c(NO3-) = -0.217 eV, from pKa(HNO3) = -1.3), but adsorbed NO3 is assigned the same correction, so it cancels in this step. A scheme that corrects free NO3- only would get dG_ads - c(NO3-).
- dE_ads_elec: the same with electronic terms only (drops the NO3- translational/rotational free energy).
- dN_e = N(*NO3) - N(*) - N(NO3-): electrons drawn from the electrode on adsorption.
- Flags: `label` is the visual inspection label of the lowest configuration (intact / physisorbed / desorbed / broken); `fin` says whether its geometry optimization finished. Where the lowest configuration is not intact and finished, `intact+fin` gives the lowest one that is.

## U = +0.0 V vs SHE

| surface | dG_ads | dE_ads_elec | dN_e | lowest configuration | label | fin | intact+fin dG_ads |
|---|---:|---:|---:|---|---|---|---:|
| ScN4-i | -0.113 | -0.621 | -0.16 | in_bi_tilted_gau_resub | intact | yes | same |
| TiN4-i | -0.904 | -1.412 | -0.42 | in_bi_metal_ecl | intact | yes | same |
| VN4-i | -1.250 | -1.759 | -0.28 | in_noo_hollow_resub | intact | yes | same |
| CrN4-i | -0.331 | -0.839 | -0.59 | in_no_bridge | intact | yes | same |
| MnN4-i | +0.062 | -0.447 | -0.46 | in_ncoord_bent | intact | yes | same |
| FeN4-i | +0.193 | -0.316 | -0.49 | in_bi_metal_atop_cc_bridge | intact | yes | same |
| CoN4-i | +0.256 | -0.252 | -0.59 | ASE_pert_from_in_bi_metal_gau | intact | yes | same |
| NiN4-i | +0.511 | +0.003 | -0.08 | in_no_bridge | physisorbed | yes | +0.662 |
| CuN4-i | +0.501 | -0.007 | -0.13 | in_bi_tilted_ecl | physisorbed | no | +0.558 |
| ZnN4-i | +0.140 | -0.368 | -0.28 | in_bi_tilted_ecl | intact | yes | same |
| YN4-i | +0.175 | -0.333 | -0.17 | in_ooo_hollow_bent | intact | yes | same |
| ZrN4-i | -0.531 | -1.039 | -0.20 | in_bi_tilted_gau | intact | yes | same |
| NbN4-i | -1.054 | -1.562 | -0.26 | in_bi_metal_gau | intact | yes | same |
| MoN4-i | -1.746 | -2.254 | -0.31 | in_noo_hollow_resub | intact | yes | same |
| TcN4-i | -1.137 | -1.646 | -0.91 | in_bi_tilted_ecl | intact | no | -1.137 |
| RuN4-i | -0.325 | -0.833 | -0.91 | in_bi_tilted_gau_resub | intact | yes | same |
| RhN4-i | +0.421 | -0.087 | -0.90 | in_bi_metal_gau | intact | yes | same |
| PdN4-i | +0.508 | -0.000 | -0.08 | in_no_bridge | physisorbed | yes | n/a |
| AgN4-i | +0.473 | -0.035 | -0.10 | in_ooo_hollow_bent | desorbed | no | n/a |
| CdN4-i | -0.111 | -0.619 | -0.22 | ASE_pert_from_in_mono_metal_atop | intact | yes | same |
| HfN4-i | -0.863 | -1.371 | -0.20 | in_bi_tilted_gau | intact | yes | same |
| TaN4-i | -1.526 | -2.034 | -0.20 | in_bi_metal_ecl | intact | yes | same |
| WN4-i | -2.248 | -2.757 | -0.71 | in_bi_metal_gau | intact | yes | same |
| ReN4-i | -1.699 | -2.207 | -0.97 | in_bi_metal_ecl | intact | yes | same |
| OsN4-i | -0.559 | -1.067 | -1.20 | in_bi_tilted_gau | intact | yes | same |
| IrN4-i | +0.507 | -0.001 | -0.03 | in_noo_hollow | physisorbed | yes | +0.591 |
| PtN4-i | +0.507 | -0.002 | -0.08 | in_no_bridge | physisorbed | yes | n/a |
| AuN4-i | +0.495 | -0.013 | -0.10 | in_noo_hollow | physisorbed | no | n/a |
| HgN4-i | +0.512 | +0.004 | +0.00 | in_mono_metal_ecl | physisorbed | yes | +0.515 |

## U = -1.0 V vs SHE

| surface | dG_ads | dE_ads_elec | dN_e | lowest configuration | label | fin | intact+fin dG_ads |
|---|---:|---:|---:|---|---|---|---:|
| ScN4-i | +0.001 | -0.508 | -0.04 | in_bi_metal_ecl | intact | yes | same |
| TiN4-i | -1.093 | -1.601 | +1.06 | in_ncoord_bent | intact | yes | same |
| VN4-i | -0.988 | -1.496 | -0.34 | in_bi_metal_gau | intact | yes | same |
| CrN4-i | +0.190 | -0.318 | -0.59 | in_bi_metal_gau | intact | yes | same |
| MnN4-i | +0.407 | -0.102 | -0.01 | in_ncoord_bent | intact | yes | same |
| FeN4-i | +0.393 | -0.115 | -0.05 | in_bi_metal_ecl | intact | yes | same |
| CoN4-i | +0.491 | -0.018 | +0.04 | in_no_bridge | physisorbed | yes | +0.658 |
| NiN4-i | +0.517 | +0.009 | +0.11 | in_mono_metal_ecl | physisorbed | yes | n/a |
| CuN4-i | +0.519 | +0.011 | +0.09 | in_mono_metal_gau | physisorbed | yes | +0.632 |
| ZnN4-i | +0.354 | -0.155 | -0.08 | in_bi_metal_ecl | intact | yes | same |
| YN4-i | +0.310 | -0.198 | -0.16 | in_bi_tilted_gau | intact | yes | same |
| ZrN4-i | -1.137 | -1.645 | +1.44 | in_ooo_hollow_bent | intact | yes | same |
| NbN4-i | -1.556 | -2.064 | +0.81 | in_ooo_hollow_bent | intact | no | -1.555 |
| MoN4-i | -1.440 | -1.949 | -0.08 | in_noo_hollow | intact | yes | same |
| TcN4-i | -0.444 | -0.953 | -0.65 | in_mono_metal_atop | intact | yes | same |
| RuN4-i | +0.382 | -0.127 | -0.99 | in_mono_metal_atop | intact | yes | same |
| RhN4-i | +0.483 | -0.025 | -0.09 | in_bi_tilted_gau | physisorbed | yes | n/a |
| PdN4-i | +0.510 | +0.001 | +0.09 | in_bi_tilted_ecl | physisorbed | no | n/a |
| AgN4-i | +0.519 | +0.010 | +0.09 | in_mono_metal_ecl | physisorbed | yes | n/a |
| CdN4-i | +0.051 | -0.457 | -0.17 | in_bi_tilted_gau | intact | yes | same |
| HfN4-i | -1.650 | -2.159 | +1.42 | in_ooo_hollow_bent | intact | no | -1.650 |
| TaN4-i | -2.837 | -3.345 | +0.02 | in_ooo_hollow_bent | intact | yes | same |
| WN4-i | -1.973 | -2.482 | +0.04 | in_ooo_hollow_bent | intact | yes | same |
| ReN4-i | -0.938 | -1.446 | -0.72 | in_bi_metal_gau | intact | yes | same |
| OsN4-i | +0.504 | -0.004 | -0.12 | in_bi_tilted_ecl | physisorbed | no | +0.512 |
| IrN4-i | +0.497 | -0.011 | -0.12 | in_mono_metal_ecl | physisorbed | yes | n/a |
| PtN4-i | +0.505 | -0.003 | +0.10 | in_mono_metal_ecl | physisorbed | yes | n/a |
| AuN4-i | +0.553 | +0.044 | +0.06 | in_mono_metal_ecl | physisorbed | yes | n/a |
| HgN4-i | +0.566 | +0.057 | -0.22 | in_mono_metal_atop | physisorbed | yes | n/a |
