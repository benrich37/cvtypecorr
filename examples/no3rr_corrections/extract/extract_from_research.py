"""Copy the inputs of this example out of the NO3_ads_gen project and record the project's own results.

This is the only file in the example that needs the `research` repo (~/vs/research) and the external drive. It is
kept here as the provenance of everything in ../data/; the example itself runs without it.

What it writes (all energies in eV):
  data/settings.json            constants and conventions the energies depend on, plus where the data came from
  data/molecules.json           the free molecules used by the fit and the steps
  data/structures.json          clean slabs (all surfaces), every cached *NO3 configuration (all surfaces) and every
                                configuration of the other adsorbates in the example steps on EXAMPLE_SURF
  data/research_reference.json  corrections, step energies and NO3 adsorption energies computed by the research repo's
                                own code, which the example scripts must reproduce

Each record keeps the entry of the energy cache (NO3_ads_gen/anl/mnc/data.json) verbatim under "cache_entry", with the
key path it came from ("cache_key") and the calculation directory relative to the project root ("calc_dir"), so any
number can be traced to PBED30_GBRV_20Ha/<calc_dir>/ion_opt/jdftx_run/out on the drive or on Anvil.

Run from anywhere, in the research conda env, with the drive mounted:
    python examples/no3rr_corrections/extract/extract_from_research.py
"""
import contextlib
import hashlib
import io
import json
import os
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from cvtypecorr.fit.rigid import Collection
from research.projects.no3_ads_gen import project_root
from research.projects.no3_ads_gen.corrections import get_collection
from research.projects.no3_ads_gen.energy_caching import CachedNrgAnl, kb
import research.projects.no3_ads_gen.no3_pathway_sim_tools as pst

RESEARCH = Path("~/vs/research").expanduser()
CACHE_DIR = RESEARCH / "NO3_ads_gen/anl/mnc"
STEP_MAP_FILE = RESEARCH / "NO3_ads_gen/anl/mnc/no3rr_pathway_map.py"
LABEL_FILE = Path(os.environ.get("HT_DATA_ANL", "~/ht_data_anl")).expanduser() / "ht_labels.json"
ANVIL_ROOT = "anvil:/anvil/projects/x-chm260113/Data/NO3_ads_gen/PBED30_GBRV_20Ha"
OUT = Path(__file__).resolve().parents[1] / "data"

# Settings of the current pathway analysis (NO3_ads_gen/anl/mnc/no3rr_pathway_map.py)
PH = 6.3
EXPECTED_TAGS = ["elec", "trans", "rot"]
ENFORCED_TAGS = ["elec"]
BIASES = ["0V", "n1p0V"]

EXAMPLE_SURF = "FeN4-i"
EXAMPLE_ADS = ["NO3", "NO3H", "O", "OH", "NO2", "H", "NO2H", "NO", "NOH", "NHO", "NHOH", "N", "NH", "NH2", "NH3"]

# cvtypecorr species name -> molecule name in the cache (mols/<name>/ in the calc tree). "X(g)" is the gas-phase twin
# of a species the chain corrects through a gas-phase reaction.
MOLECULES = {
    "H2": "H2_gas", "H2O": "H2O", "H3O+": "H3O_p1", "OH-": "OH_n1",
    "NO3H": "NO3H", "NO3-": "NO3_n1", "NO2H": "NO2H", "NO2-": "NO2_n1",
    "NO": "NO_gas", "NH3": "NH3", "NH4+": "NH4_p1",
    "NHO": "NHO", "NH2O": "NH2O", "NH2OH": "NH2OH",
    **{f"{s}(g)": f"{s}_gas" for s in ["NH3", "NH2", "NH", "N", "NHO", "NOH", "NH2O", "NHOH", "NH2OH", "NH3O"]},
}

# Example steps: (initial state, final state) keys of intermediate_next_step_map in no3rr_pathway_map.py
EXAMPLE_STEPS = [
    ("*", "NO3"), ("*", "O"), ("*", "NO2"), ("*", "H"),
    ("NO3", "NO3H"), ("NO3", "O"), ("NO3H", "NO2"),
    ("O", "OH"), ("OH", "*"), ("H", "*"),
    ("NO2", "NO2H"), ("NO2H", "NO"), ("NO", "NHO"), ("NO", "NOH"), ("NHO", "NHOH"), ("NHOH", "NH"),
    ("N", "NH"), ("NH", "NH2"), ("NH2", "NH3"),
]


def git_head(repo: Path) -> str:
    return subprocess.run(["git", "-C", str(repo), "rev-parse", "--short", "HEAD"],
                          capture_output=True, text=True, check=True).stdout.strip()


def git_is_clean(repo: Path, paths: list[str]) -> bool:
    out = subprocess.run(["git", "-C", str(repo), "status", "--porcelain", "--", *paths],
                         capture_output=True, text=True, check=True).stdout
    return not out.strip()


def sha1(path: Path) -> str:
    return hashlib.sha1(path.read_bytes()).hexdigest()


def file_info(path: Path) -> dict | None:
    if not path.exists():
        return None
    st = path.stat()
    return {"size": st.st_size,
            "mtime": datetime.fromtimestamp(st.st_mtime, timezone.utc).isoformat(timespec="seconds")}


def finite_or_none(obj):
    """NaN/inf -> None (JSON null), recursively; numpy scalars -> float."""
    if isinstance(obj, dict):
        return {k: finite_or_none(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [finite_or_none(v) for v in obj]
    if isinstance(obj, (float, np.floating)):
        return float(obj) if np.isfinite(obj) else None
    return obj


def make_cna() -> CachedNrgAnl:
    cna = CachedNrgAnl(expected_tags=EXPECTED_TAGS, enforced_tags=ENFORCED_TAGS, save_dir=CACHE_DIR)
    with contextlib.redirect_stdout(io.StringIO()):
        cna.corrs = get_collection(cna=cna, pH=PH)
    return cna


def structure_record(cna, labels, surf, ads, bias, config) -> dict:
    """One configuration (ads == "*" is the clean slab), with its cache entry and where it lives in the calc tree."""
    if ads == "*":
        cache_key = ["structures", surf, "surfs", bias, "clean"]
        calc_dir = f"surfs/{surf}/{bias}"
    else:
        cache_key = ["structures", surf, ads, bias, config]
        calc_dir = f"ads/{surf}/{ads}/{bias}/{config}"
    root = project_root / calc_dir
    out_file = root / "ion_opt/jdftx_run/out"
    label = labels.get(str(out_file), {})
    entry = cna.data
    for k in cache_key:
        entry = entry[k]
    return {
        "surf": surf, "ads": ads, "bias": bias, "config": config,
        "calc_dir": calc_dir,
        "label": {k: label[k] for k in ("state", "geometry", "pseudonyms", "human_approved") if k in label},
        # finished.txt is the project's marker for a converged ion_opt; calcs run by Sophie carry their own marker
        "ion_opt_finished": (root / "ion_opt/finished.txt").exists() or (root / "sophie_calc_finished.txt").exists(),
        "from_sophie": (root / "sophie_mark.txt").exists(),
        "out_file": file_info(out_file),
        "cache_key": cache_key,
        "cache_entry": entry,
    }


def molecule_record(cna, species, name) -> dict:
    root = project_root / "mols" / name
    return {
        "species": species, "name": name,
        "calc_dir": f"mols/{name}",
        "out_file": file_info(root / "ion_opt/jdftx_run/out"),
        "vib_out_file": file_info(root / "vib/jdftx_run/out"),
        "cache_key": ["molecules", name],
        "cache_entry": cna.data["molecules"][name],
    }


def load_step_funcs() -> dict:
    """The "func" lambda source of each example step, read from the step map so the reference uses it verbatim."""
    src = STEP_MAP_FILE.read_text()
    body = src[src.index("intermediate_next_step_map = {"):]
    funcs, cur, sub = {}, None, None
    for line in body.splitlines()[1:]:
        if line.startswith("}"):
            break
        if m := re.match(r'^    "([^"]+)": \{', line):
            cur = m.group(1)
        elif m := re.match(r'^        "([^"]+)": \{', line):
            sub = m.group(1)
        elif m := re.match(r'^\s+"func": (lambda .*),\s*$', line):
            funcs[(cur, sub)] = m.group(1)
    return funcs


def research_steps(cna, funcs) -> list[dict]:
    ns = {**vars(pst), "pH": PH}
    zero = Collection(corrections={k: 0.0 for k in cna.corrs.corrections})
    real = cna.corrs
    rows = []
    for key in EXAMPLE_STEPS:
        func = eval(funcs[key], ns)
        for bias in BIASES:
            cna.corrs = real
            g_corr, s_a, s_b = func(EXAMPLE_SURF, cna, bias)
            cna.corrs = zero
            g_raw, _, _ = func(EXAMPLE_SURF, cna, bias)
            rows.append({"step": f"{key[0]} -> {key[1]}", "bias": bias, "dG_raw": g_raw, "dG_corrected": g_corr,
                         "configs_returned": [s_a, s_b], "step_map_func": funcs[key]})
    cna.corrs = real
    return rows


def research_no3_adsorption(cna, surfs) -> list[dict]:
    rows = []
    for surf in surfs:
        for bias in BIASES:
            g, _, config = pst.get_generic_no3_ads_nrg(cna, surf, bias, "*", "NO3", rct_tuples=[], pdt_tuples=[])
            rows.append({"surf": surf, "bias": bias, "dG": g, "config": config})
    return rows


def main():
    for p in (project_root, LABEL_FILE, STEP_MAP_FILE):
        assert p.exists(), f"missing {p} (is the drive mounted?)"
    code_files = ["src/research/projects/no3_ads_gen/corrections.py",
                  "src/research/projects/no3_ads_gen/energy_caching.py",
                  "src/research/projects/no3_ads_gen/no3_pathway_sim_tools.py",
                  "NO3_ads_gen/anl/mnc/no3rr_pathway_map.py",
                  "NO3_ads_gen/anl/mnc/data.json"]
    pst.print_errors = False
    cna = make_cna()
    labels = json.loads(LABEL_FILE.read_text())
    surfs = sorted(cna.data["structures"])
    OUT.mkdir(parents=True, exist_ok=True)

    settings = {
        "_about": "Constants and conventions of the NO3_ads_gen energies as used by the current pathway analysis "
                  "(no3rr_pathway_map.py). Energies in eV.",
        "T_K": cna.T,
        "kB_eV_per_K": kb,
        "v0_SHE_eV": cna.v0,
        "pH": PH,
        "energy_tags": EXPECTED_TAGS,
        "required_tags": ENFORCED_TAGS,
        "mu_eV": {"No_bias": 0.0, **{b: cna.data["mu_dict"][b] for b in BIASES}},
        "bias_vs_SHE_V": {b: -cna.data["mu_dict"][b] - cna.v0 for b in BIASES},
        "conc_M": {"fit": {"H2O": 55.5}, "steps": {"H2O": 55.34, "H3O+": 10 ** -PH}},
        "provenance": {
            "extracted": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "research_commit": git_head(RESEARCH),
            "research_files_unmodified_at_commit": git_is_clean(RESEARCH, code_files),
            "research_file_sha1": {f: sha1(RESEARCH / f) for f in code_files},
            "energy_cache": "research:NO3_ads_gen/anl/mnc/data.json",
            "energy_cache_sha1": sha1(CACHE_DIR / "data.json"),
            "labels": str(LABEL_FILE).replace(str(Path.home()), "~"),
            "labels_sha1": sha1(LABEL_FILE),
            "project_root_local": str(project_root),
            "project_root_anvil": ANVIL_ROOT,
            "calc_dir_is_relative_to": "project root (PBED30_GBRV_20Ha)",
        },
    }

    molecules = {"_about": "Free molecules. cache_entry is verbatim data.json['molecules'][name]; the calculation is "
                           "<project root>/<calc_dir>/ion_opt (geometry, elec) and .../vib (vib, trans, rot terms).",
                 "records": [molecule_record(cna, s, n) for s, n in MOLECULES.items()]}

    records = [structure_record(cna, labels, surf, "*", b, "clean") for surf in surfs for b in BIASES]
    for surf in surfs:
        for ads in (EXAMPLE_ADS if surf == EXAMPLE_SURF else ["NO3"]):
            for b in BIASES:
                for config in cna.data["structures"][surf].get(ads, {}).get(b, {}):
                    records.append(structure_record(cna, labels, surf, ads, b, config))
    structures = {"_about": "Slab calculations. ads '*' is the clean slab. cache_entry is verbatim "
                            "data.json[cache_key]; the calculation is <project root>/<calc_dir>/ion_opt/jdftx_run/out. "
                            "label is the ht_browser label of that out file.",
                  "records": records}

    reference = {
        "_about": "Values computed by the research repo's own code (corrections.get_collection and the step functions "
                  "of no3_pathway_sim_tools via the step map of no3rr_pathway_map.py); the example must reproduce them.",
        "corrections": cna.corrs.corrections,
        "steps": {"surf": EXAMPLE_SURF, "rows": research_steps(cna, load_step_funcs())},
        "no3_adsorption": research_no3_adsorption(cna, surfs),
    }

    for name, obj in [("settings", settings), ("molecules", molecules), ("structures", structures),
                      ("research_reference", reference)]:
        (OUT / f"{name}.json").write_text(json.dumps(finite_or_none(obj), indent=1, allow_nan=False) + "\n")
    print(f"wrote {OUT}: {len(molecules['records'])} molecules, {len(records)} structures "
          f"({sum(r['ads'] == 'NO3' for r in records)} *NO3 configurations on {len(surfs)} surfaces)")


if __name__ == "__main__":
    main()
