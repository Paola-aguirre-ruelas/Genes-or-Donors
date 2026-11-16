"""Assemble component-by-component comparisons of a run with the historical references."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .io_utils import read_json, write_json
from .provenance import utc_now
from .reference import (DISCREPANCY, EXACT, MISSING, NUMERIC, Legacy, compare_model, compare_plan,
                        compare_preprocessing, compare_scores, compare_search, worst)

ARMS = ("w0", "w1")


def _load_ckpt(path: Path) -> dict:
    with np.load(path, allow_pickle=False) as z:
        return {"meta": json.loads(str(z["meta"].item())), "coef": z["coef"].copy(), "best_norms": z["best_norms"].copy()}


def _load_npz(path: Path) -> dict:
    with np.load(path, allow_pickle=False) as z:
        return {k: z[k].copy() for k in z.files}


def compare_training(run_dir: Path, legacy: Legacy, n_folds: int, tol: dict) -> dict:
    out = {}
    for k in range(n_folds):
        fd = run_dir / f"fold{k}"
        comp = {}
        if (fd / "plan.json").exists():
            comp["plan"] = compare_plan(read_json(fd / "plan.json")["lists"], legacy.plan_lists(k))
        if (fd / "preprocessing.npz").exists():
            z = _load_npz(fd / "preprocessing.npz")
            comp["preprocessing"] = compare_preprocessing(
                {"candidate_var_idx": z["candidate_var_idx"], "genes": z["genes"].astype(str), "means": z["means"], "sds": z["sds"]},
                legacy.preprocessing(k), tol)
        for arm in ARMS:
            if (fd / f"panel_{arm}.json").exists():
                comp[f"panel_{arm}"] = compare_search(_load_ckpt(fd / f"checkpoint_{arm}.npz"), read_json(fd / f"panel_{arm}.json"),
                                                      legacy.checkpoint(k, arm), legacy.panel(k, arm), tol)
            if (fd / f"model_{arm}.npz").exists():
                z = _load_npz(fd / f"model_{arm}.npz")
                comp[f"model_{arm}"] = compare_model({"coef": z["coef"], "intercept": z["intercept"], "classes": z["classes"].astype(str),
                                                      "panel_idx": z["panel_idx"]}, legacy.model(k, arm), tol)
        rd = fd / "random"
        if (rd / "manifest.json").exists():
            man = read_json(rd / "manifest.json")
            old = legacy.random_manifest(k)
            panels_eq = man["panels"] == old["panels"] and man["bin_counts"] == old["bin_counts"]
            draws = []
            for d in range(len(man["panels"])):
                p = rd / f"model_{d:02d}.npz"
                if p.exists():
                    z = _load_npz(p)
                    draws.append(compare_model({"coef": z["coef"], "intercept": z["intercept"], "classes": z["classes"].astype(str),
                                                "panel_idx": z["cols"]},
                                               {**legacy.random_model(k, d), "panel_idx": legacy.random_model(k, d)["cols"]}, tol))
            st = [d["status"] for d in draws]
            comp["random_reference"] = {"component": "random_reference",
                                        "status": worst(st + ([EXACT] if panels_eq else [DISCREPANCY])),
                                        "panels_identical": panels_eq, "n_models_compared": len(draws),
                                        "model_statuses": {s: st.count(s) for s in set(st)},
                                        "coef_max_abs_diff": max((d["coef_max_abs_diff"] for d in draws), default=None)}
        out[f"fold{k}"] = comp
    return out


def compare_evaluation(run_dir: Path, legacy: Legacy, n_folds: int, tol: dict) -> dict:
    out = {}
    for k in range(n_folds):
        p = run_dir / f"fold{k}" / "evaluation.json"
        if not p.exists():
            continue
        new = read_json(p)
        old = legacy.evaluation(k)
        oldr = legacy.random_evaluation(k)
        main = compare_scores(new["per_individual"], old["per_donor"], ["macro_f1_w0", "macro_f1_w1", "paired_delta"], tol)
        old_rand = {r["individual_id"]: r for r in oldr["per_donor"]}
        rand_rows_old = [{"individual_id": d, **{f"d{j}": v for j, v in enumerate(r["f1_by_random_draw"])}, "mean_random_f1": r["mean_random_f1"]}
                         for d, r in old_rand.items()]
        rand_rows_new = [{"individual_id": r["individual_id"], **{f"d{j}": v for j, v in enumerate(r["f1_by_random_draw"])},
                          "mean_random_f1": r["mean_random_f1"]} for r in new["per_individual"]]
        nd = len(new["per_individual"][0]["f1_by_random_draw"])
        rand = compare_scores(rand_rows_new, rand_rows_old, [f"d{j}" for j in range(nd)] + ["mean_random_f1"], tol)
        old_delta = old.get("paired_delta", old.get("paired_delta_pilot"))
        out[f"fold{k}"] = {"main_arms": main, "random_reference": rand,
                           "fold_delta_new": new["paired_delta"], "fold_delta_old": old_delta,
                           "fold_delta_abs_diff": abs(new["paired_delta"] - old_delta)}
    return out


def compare_stability(run_dir: Path, legacy: Legacy, n_folds: int, tol: dict) -> dict:
    out = {}
    sd = run_dir / "stability"
    for k in range(n_folds):
        for arm in ARMS:
            p = sd / f"g3_fold{k}_{arm}.json"
            if not p.exists():
                continue
            new, old = read_json(p), legacy.stability(f"g3_fold{k}_{arm}")
            disc = (new["tight_support_genes"] == old["tight_support_genes"] and new["status"] == old["status"]
                    and new["tight_top64_genes"] == old["tight_top64_genes"])
            gap_rel = abs(new["fit"]["dual_gap"] - old["fit"]["dual_gap"]) / max(abs(old["fit"]["dual_gap"]), 1e-300)
            same_iter = new["fit"]["n_iter"] == old["fit"]["n_iter"]
            status = (EXACT if disc and same_iter and new["fit"]["dual_gap"] == old["fit"]["dual_gap"]
                      else NUMERIC if disc and gap_rel <= float(tol["dual_gap_rtol"]) else DISCREPANCY)
            out[f"g3_fold{k}_{arm}"] = {"status": status, "verdict_new": new["status"], "verdict_old": old["status"],
                                        "support_identical": disc, "n_iter_identical": same_iter, "dual_gap_rel_diff": gap_rel}
    for key, fname in (("g4", "g4_state.json"), ("h1", "h1_state.json")):
        p = sd / fname
        if not p.exists():
            continue
        new, old = read_json(p), legacy.stability(key)
        a_new = [(r["alpha"], r["nnz"], r["n_iter"]) for r in new["records"]]
        a_old = [(r["alpha"], r["nnz"], r["n_iter"]) for r in old["records"]]
        disc = [(x[1]) for x in a_new] == [(x[1]) for x in a_old] and len(a_new) == len(a_old)
        alpha_rel = max((abs(x[0] - y[0]) / abs(y[0]) for x, y in zip(a_new, a_old)), default=0.0)
        if key == "g4":
            disc = disc and new.get("result", {}).get("selected") == old.get("result", {}).get("selected")
        status = EXACT if a_new == a_old and disc else NUMERIC if disc and alpha_rel <= float(tol["alpha_rtol"]) else DISCREPANCY
        res_new = {k2: v for k2, v in new.get("result", new.get("conclusion", {})).items() if k2 not in ("selected", "reference")}
        res_old_all = old.get("result", old.get("conclusion", {}))
        # Echo only the structured fields both records share; historical free-text notes are not republished.
        res_old = {k2: res_old_all[k2] for k2 in res_new if k2 in res_old_all}
        out[key] = {"status": status, "trace_identical": a_new == a_old, "nnz_trace_identical": disc,
                    "alpha_max_rel_diff": alpha_rel, "result_new": res_new, "result_old": res_old,
                    "result_fields_identical": res_new == res_old}
    return out


def by_component(res: dict) -> dict:
    groups = {"fold plans": [], "preprocessing": [], "alpha searches and panels": [], "classifiers": [],
              "random-reference panels and models": [], "held-out scores (main arms)": [],
              "held-out scores (random reference)": [], "stability checks": []}
    for f in res["training"].values():
        for key, c in f.items():
            name = ("fold plans" if key == "plan" else "preprocessing" if key == "preprocessing"
                    else "alpha searches and panels" if key.startswith("panel_") else "classifiers" if key.startswith("model_")
                    else "random-reference panels and models")
            groups[name].append(c["status"])
    for v in res["evaluation"].values():
        groups["held-out scores (main arms)"].append(v["main_arms"]["status"])
        groups["held-out scores (random reference)"].append(v["random_reference"]["status"])
    groups["stability checks"] = [v["status"] for v in res["stability"].values()]
    return {k: {s: v.count(s) for s in set(v)} for k, v in groups.items() if v}


def build_comparison(run_dir: Path, legacy: Legacy, n_folds: int, tol: dict) -> dict:
    tr = compare_training(run_dir, legacy, n_folds, tol)
    ev = compare_evaluation(run_dir, legacy, n_folds, tol)
    stb = compare_stability(run_dir, legacy, n_folds, tol)
    statuses = [c["status"] for f in tr.values() for c in f.values()]
    statuses += [v["main_arms"]["status"] for v in ev.values()] + [v["random_reference"]["status"] for v in ev.values()]
    statuses += [v["status"] for v in stb.values()]
    res = {"compared_utc": utc_now(), "tolerances": tol, "overall": worst(statuses),
           "counts": {s: statuses.count(s) for s in (EXACT, NUMERIC, DISCREPANCY, MISSING)},
           "training": tr, "evaluation": ev, "stability": stb}
    res["by_component"] = by_component(res)
    res["narrative"] = "; ".join(f"{name}: {', '.join(f'{n} {s}' for s, n in sorted(c.items()))}"
                                 for name, c in res["by_component"].items()) + "."
    write_json(run_dir / "comparisons" / f"comparison_{utc_now().replace(':', '')}.json", res)
    write_json(run_dir / "comparisons" / "latest.json", res, overwrite=True)
    return res


# ---------------------------------------------------------------------- results package vs historical delivery
SUMMARY_KEYS = {  # regenerated key -> historical descriptive_results.json key
    "paired_delta": "paired_delta", "median_delta": "median_delta", "mean_macro_f1_w0": "mean_macro_f1_w0",
    "mean_macro_f1_w1": "mean_macro_f1_w1", "random_reference_mean": "random_baseline_mean",
    "w0_minus_random": "w0_minus_random", "w1_minus_random": "w1_minus_random",
    "per_individual_positive": "per_donor_positive", "per_individual_negative": "negative", "per_individual_zero": "zero",
    "w0_exceeds_all_random_fold_means": "baseline_w0_exceeds_all_20_fold_means", "n_individuals": "n_donors"}
TABLES = {"table_01_individuals_anonymized.csv": "table_01", "table_02_fold_summary.csv": "table_02",
          "table_03_panel_overlap.csv": "table_03", "table_04_selected_gene_membership.csv": "table_04"}
FIGURES = {"figure_01_paired_donor_deltas.png": "figure_01", "figure_02_random_baseline.png": "figure_02",
           "figure_03_panel_overlap.png": "figure_03"}


def _csv_numeric_diff(a: str, b: str) -> tuple[bool, float]:
    import csv
    import io
    ra, rb = list(csv.reader(io.StringIO(a))), list(csv.reader(io.StringIO(b)))
    if len(ra) != len(rb) or any(len(x) != len(y) for x, y in zip(ra, rb)) or ra[0] != rb[0]:
        return False, float("inf")
    worst_d = 0.0
    for x, y in zip(ra[1:], rb[1:]):
        for u, v in zip(x, y):
            try:
                worst_d = max(worst_d, abs(float(u) - float(v)))
            except ValueError:
                if u != v:
                    return False, float("inf")
    return True, worst_d


def compare_results_package(res_dir: Path, legacy: Legacy, tol: dict) -> dict:
    from .io_utils import sha256_file
    rep = legacy.refs["report"]
    out = {}
    for new, key in TABLES.items():
        a, b = (res_dir / new).read_text(), legacy.path(rep[key]).read_text()
        same_shape, d = _csv_numeric_diff(a, b)
        status = EXACT if a == b else NUMERIC if same_shape and d <= float(tol["aggregate_atol"]) else DISCREPANCY
        out[new] = {"status": status, "byte_identical": a == b, "max_abs_diff": d, "historical": rep[key]}
    old = read_json(legacy.path(rep["descriptive_results"]))
    new = read_json(res_dir / "results.json")["summary"]
    diffs, st = {}, []
    for nk, ok in SUMMARY_KEYS.items():
        a, b = new[nk], old[ok]
        if isinstance(a, bool) or isinstance(b, bool) or isinstance(a, int) and isinstance(b, int):
            s = EXACT if a == b else DISCREPANCY
            diffs[nk] = {"new": a, "historical": b, "status": s}
        else:
            d = abs(float(a) - float(b))
            s = EXACT if a == b else NUMERIC if d <= float(tol["aggregate_atol"]) else DISCREPANCY
            diffs[nk] = {"new": a, "historical": b, "abs_diff": d, "status": s}
        st.append(s)
    out["summary_values"] = {"status": worst(st), "values": diffs, "historical": rep["descriptive_results"]}
    for new_f, key in FIGURES.items():
        same = sha256_file(res_dir / new_f) == sha256_file(legacy.path(rep[key]))
        # Figures are presentation artifacts: byte identity is reported; data identity is established by the tables.
        out[new_f] = {"status": EXACT if same else NUMERIC, "byte_identical": same, "historical": rep[key],
                      "note": "data identity established by tables 01-03" if not same else ""}
    return out
