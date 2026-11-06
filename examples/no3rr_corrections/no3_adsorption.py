"""Step 3: NO3- adsorption free energies on every MN4-i surface, at 0 and -1 V vs SHE (the comparison set).

    dG_ads = G(*NO3) - G(*) - G(NO3-)          step "* -> NO3" of the pathway analysis
    correction = c(NO3-) [*NO3 carries the NO3- label] - c(NO3-) [free] = 0, so dG_ads is also the corrected value.

G(*NO3) is the lowest of all cached *NO3 configurations (as in the pathway analysis); dG_ads_intact_finished repeats
the minimum over configurations labelled 'intact' whose ion_opt finished. dE_ads_elec uses the "elec" term only (no
trans/rot free energy of NO3-). dN_e = N(*NO3) - N(*) - N(NO3-) is the number of electrons drawn from the electrode.
ion_opt_finished counts ion_opt/finished.txt, or sophie_calc_finished.txt for calcs marked with sophie_mark.txt (run by
Sophie; column sophie_mark).

Run: python no3_adsorption.py
  -> output/no3_adsorption_energies.csv     one row per surface and bias
     output/no3_adsorption_all_configs.csv  one row per configuration
     output/no3_adsorption_energies.md      the per-surface table with definitions, readable on its own
"""
import csv
import math

from pymatgen.core import Element

from energies import HERE, MOLECULES, MU, REFERENCE, SETTINGS, STRUCTURES, configs, grand_free_energy, lowest
from fit_corrections import fit_corrections

BIASES = ["0V", "n1p0V"]
ELEC = ("elec",)


def periodic_order(surf: str) -> tuple[int, int]:
    metal = Element(surf.split("N4")[0])
    return metal.row, metal.group


def main():
    corrections = fit_corrections()[0].corrections
    corr = corrections["NO3-"] - corrections["NO3-"]  # *NO3 is labelled NO3-: the correction cancels
    no3 = MOLECULES["NO3-"]["cache_entry"]
    surfs = sorted({r["surf"] for r in STRUCTURES}, key=periodic_order)
    research = {(r["surf"], r["bias"]): r for r in REFERENCE["no3_adsorption"]}

    per_surface, per_config, max_err = [], [], 0.0
    for surf in surfs:
        for bias in BIASES:
            mu = MU[bias]
            G_clean, clean = lowest(surf, "*", bias)
            G_no3 = grand_free_energy(no3, mu)
            dE_ref = grand_free_energy(clean["cache_entry"], mu, ELEC) + grand_free_energy(no3, mu, ELEC)
            rows = []
            for rec in configs(surf, "NO3", bias):
                G = grand_free_energy(rec["cache_entry"], mu)
                rows.append({
                    "surface": surf, "bias": bias, "config": rec["config"],
                    "dG_ads_eV": G - G_clean - G_no3 + corr,
                    "dE_ads_elec_eV": grand_free_energy(rec["cache_entry"], mu, ELEC) - dE_ref,
                    "dN_e": rec["cache_entry"]["N"] - clean["cache_entry"]["N"] - no3["N"],
                    "G_slab_NO3_eV": G,
                    "label_state": rec["label"].get("state", ""), "label_geometry": rec["label"].get("geometry", ""),
                    "ion_opt_finished": rec["ion_opt_finished"], "sophie_mark": rec["from_sophie"],
                    "calc_dir": rec["calc_dir"],
                    "out_size": (rec["out_file"] or {}).get("size"), "out_mtime": (rec["out_file"] or {}).get("mtime"),
                })
            usable = [r for r in rows if math.isfinite(r["dG_ads_eV"])]  # drop NaN (null terms in the cache)
            best = min(usable, key=lambda r: r["dG_ads_eV"])
            good = [r for r in usable if r["label_state"] == "intact" and r["ion_opt_finished"]]
            best_good = min(good, key=lambda r: r["dG_ads_eV"]) if good else None
            for r in rows:
                r["is_lowest"] = r is best
            per_config += rows

            ref = research[(surf, bias)]
            assert ref["config"] == best["config"], (surf, bias, ref["config"], best["config"])
            max_err = max(max_err, abs(ref["dG"] - best["dG_ads_eV"]))
            per_surface.append({
                "surface": surf, "bias": bias, "U_vs_SHE_V": round(SETTINGS["bias_vs_SHE_V"][bias], 6),
                "dG_ads_eV": best["dG_ads_eV"], "dE_ads_elec_eV": best["dE_ads_elec_eV"], "dN_e": best["dN_e"],
                "config": best["config"], "label_state": best["label_state"],
                "label_geometry": best["label_geometry"], "ion_opt_finished": best["ion_opt_finished"],
                "sophie_mark": best["sophie_mark"],
                "n_configs": len(rows),
                "dG_ads_intact_finished_eV": best_good["dG_ads_eV"] if best_good else None,
                "config_intact_finished": best_good["config"] if best_good else "",
                "G_slab_NO3_eV": best["G_slab_NO3_eV"], "G_slab_eV": G_clean,
                "G_NO3-_eV": G_no3, "calc_dir": best["calc_dir"], "clean_calc_dir": clean["calc_dir"],
            })
    assert max_err < 1e-8, f"NO3 adsorption energies differ from the research repo by {max_err:.2e} eV"

    out = HERE / "output"
    for name, rows in [("no3_adsorption_energies", per_surface), ("no3_adsorption_all_configs", per_config)]:
        with open(out / f"{name}.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    write_markdown(per_surface, corrections["NO3-"], out / "no3_adsorption_energies.md")
    print(f"{len(per_surface)} surface/bias rows, {len(per_config)} configurations; "
          f"matches research get_generic_no3_ads_nrg (value and configuration): max |diff| = {max_err:.1e} eV")
    print("wrote output/no3_adsorption_energies.{csv,md}, output/no3_adsorption_all_configs.csv")


def write_markdown(rows: list[dict], c_no3: float, path):
    s = SETTINGS
    no3 = MOLECULES["NO3-"]["cache_entry"]
    g_tr = sum(no3["E"][t] - no3["TS"][t] for t in ("trans", "rot"))
    lines = [
        "# NO3- adsorption free energies on MN4-i single-atom sites",
        "",
        "Reaction `* + NO3-(aq) -> *NO3` at fixed electrode potential, from the NO3_ads_gen calculations (PBE-D3(0), "
        "GBRV, 20 Ha, CANDLE LinearPCM, 0.1 M NaF, grand canonical at fixed U). Energies in eV. Written by "
        "`no3_adsorption.py`; every number is in `no3_adsorption_energies.csv` with its components and calc dirs "
        "(`calc_dir`, relative to `PBED30_GBRV_20Ha`, on Anvil at "
        f"`{s['provenance']['project_root_anvil']}`), and every configuration in `no3_adsorption_all_configs.csv`.",
        "",
        "**Definition.** dG_ads = G(*NO3) - G(*) - G(NO3-), each G a grand free energy at the electrode "
        f"potential, G = E - TS - N*mu with mu = -(U_SHE + {s['v0_SHE_eV']}) eV, so the electrons drawn from the "
        "electrode (dN_e) are counted in G and no separate electron term appears.",
        "",
        "- G(*), G(*NO3): DFT free energy of the slab (E - TS of the electronic smearing). No vibrations or "
        "zero-point energy. G(*NO3) is the lowest of all cached *NO3 configurations.",
        f"- G(NO3-): DFT free energy of solvated NO3- plus its solvated translational and rotational free energy "
        f"({g_tr:+.3f} eV at {s['T_K']:.0f} K), at 1 M. No vibrations or zero-point energy.",
        "- Correction: none in effect. NO3- carries an empirical correction "
        f"(c(NO3-) = {c_no3:+.3f} eV, from pKa(HNO3) = -1.3), but adsorbed NO3 is "
        "assigned the same correction, so it cancels in this step. A scheme that corrects free NO3- only would get "
        "dG_ads - c(NO3-).",
        "- dE_ads_elec: the same with electronic terms only (drops the NO3- translational/rotational free energy).",
        "- dN_e = N(*NO3) - N(*) - N(NO3-): electrons drawn from the electrode on adsorption.",
        "- Flags: `label` is the visual inspection label of the lowest configuration (intact / physisorbed / desorbed "
        "/ broken); `fin` says whether its geometry optimization finished. Where the lowest configuration is not "
        "intact and finished, `intact+fin` gives the lowest one that is.",
        "",
    ]
    for bias in BIASES:
        lines += [f"## U = {s['bias_vs_SHE_V'][bias]:+.1f} V vs SHE", "",
                  "| surface | dG_ads | dE_ads_elec | dN_e | lowest configuration | label | fin | intact+fin dG_ads |",
                  "|---|---:|---:|---:|---|---|---|---:|"]
        for r in (r for r in rows if r["bias"] == bias):
            alt = r["dG_ads_intact_finished_eV"]
            same = r["config_intact_finished"] == r["config"]
            lines.append(f"| {r['surface']} | {r['dG_ads_eV']:+.3f} | {r['dE_ads_elec_eV']:+.3f} | {r['dN_e']:+.2f} | "
                         f"{r['config']} | {r['label_state']} | {'yes' if r['ion_opt_finished'] else 'no'} | "
                         f"{'same' if same else ('n/a' if alt is None else f'{alt:+.3f}')} |")
        lines.append("")
    path.write_text("\n".join(lines))


if __name__ == "__main__":
    main()
