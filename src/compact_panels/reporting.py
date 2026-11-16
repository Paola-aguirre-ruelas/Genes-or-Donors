"""Programmatic results package: tables, figures, machine-readable results, SHA-256 manifest.

Everything is derived from the run's recorded outputs (``fold<k>/evaluation.json``,
``fold<k>/panel_<arm>.json`` and ``stability/*.json``); no number is typed by hand.
Source individual identifiers are replaced by within-fold ordinal codes (``F<k>_D<nn>``,
ordered by the sorted source identifier) in every public table. These codes are pseudonyms: with
the public data the fold assignment and ordering can be recomputed (see docs/DATA.md).

Aggregation follows the historical summary code operation by operation, so that
identical inputs give bit-identical summaries.
"""
from __future__ import annotations

import csv
import io
import math
from pathlib import Path

import numpy as np

from .io_utils import read_json, require, sha256_bytes, sha256_file, write_bytes, json_bytes
from .provenance import environment, utc_now

FIGURE_NAMES = ("figure_01_paired_donor_deltas.png", "figure_02_random_baseline.png", "figure_03_panel_overlap.png")
GREY, W0C, W1C = "#8a8f98", "#1f5f99", "#c4581b"
FOLD_COLORS = ("#1f5f99", "#c4581b", "#3f8f4f", "#8e3b8e")


def _close(a, b, label, atol=1e-10):
    require(math.isclose(float(a), float(b), rel_tol=0.0, abs_tol=atol), f"{label} differs: {a} vs {b}")


def collect(run_dir: Path, cfg) -> tuple[list, list, list, dict]:
    n_folds = cfg.n_folds
    n_draws = int(cfg["random_reference"]["draws"])
    k_genes = cfg.k
    paired, byfold, panel_rows, used = [], [], [], set()
    sign = {"positive": 0, "negative": 0, "zero": 0}
    beats_all = 0
    sources = {}
    for k in range(n_folds):
        ep = run_dir / f"fold{k}" / "evaluation.json"
        require(ep.exists(), f"fold {k} has not been scored")
        sources[f"fold{k}/evaluation.json"] = sha256_file(ep)
        ev = read_json(ep)
        a = {str(r["individual_id"]): r for r in ev["per_individual"]}
        require(not (used & set(a)), f"fold {k}: individual appears in several folds")
        used.update(a)
        fold_pairs, fold_rand, fold_w0, fold_w1, draw_matrix = [], [], [], [], []
        for i, person in enumerate(sorted(a), start=1):
            x = a[person]
            w0, w1 = float(x["macro_f1_w0"]), float(x["macro_f1_w1"])
            delta, random_mean = float(x["paired_delta"]), float(x["mean_random_f1"])
            draws = np.asarray(x["f1_by_random_draw"], dtype=np.float64)
            require(draws.shape == (n_draws,) and np.isfinite(draws).all(), f"fold {k}: invalid random draws")
            _close(delta, w0 - w1, f"fold {k} paired delta")
            _close(random_mean, draws.mean(), f"fold {k} random mean")
            label = "positive" if delta > 1e-12 else ("negative" if delta < -1e-12 else "zero")
            sign[label] += 1
            beats_all += int(bool(np.all(w0 > draws)))
            paired.append(dict(donor_code=f"F{k}_D{i:02d}", fold=k, n_test_cells=int(x["n_test_cells"]),
                               n_evaluable_classes=int(x["n_evaluable_classes"]), macro_f1_w0=w0, macro_f1_w1=w1,
                               paired_delta=delta, mean_random_f1=random_mean, w0_minus_mean_random=w0 - random_mean))
            fold_pairs.append(delta); fold_rand.append(random_mean)
            fold_w0.append(w0); fold_w1.append(w1); draw_matrix.append(draws)
        draw_matrix = np.asarray(draw_matrix, dtype=float)
        draws_by_fold = draw_matrix.mean(axis=0)
        _close(np.mean(fold_pairs), ev["paired_delta"], f"fold {k} delta")
        byfold.append(dict(fold=k, donors=len(a), test_cells=int(ev["n_test_cells"]),
                           macro_f1_w0=float(np.mean(fold_w0)), macro_f1_w1=float(np.mean(fold_w1)),
                           paired_delta=float(np.mean(fold_pairs)), mean_random_f1=float(np.mean(fold_rand)),
                           random_min_draw=float(np.min(draws_by_fold)), random_max_draw=float(np.max(draws_by_fold)),
                           w0_minus_random=float(np.mean(fold_w0) - np.mean(fold_rand)), draws=draws_by_fold.tolist()))
        for arm in ("w0", "w1"):
            pp = run_dir / f"fold{k}" / f"panel_{arm}.json"
            sources[f"fold{k}/panel_{arm}.json"] = sha256_file(pp)
            obj = read_json(pp)
            genes = list(obj["panel_genes"])
            require(len(genes) == k_genes and len(set(genes)) == k_genes, f"fold {k} {arm}: invalid panel")
            panel_rows.append({"fold": k, "arm": arm, "genes": set(genes), "exact64": bool(obj["exact64"]),
                               "discarded": int(obj["discarded"]), "alpha_final": float(obj["alpha_final"]),
                               "n_fits": int(obj["n_fits"])})
    all_d = np.asarray([r["paired_delta"] for r in paired], dtype=np.float64)
    all_w0 = np.asarray([r["macro_f1_w0"] for r in paired], dtype=np.float64)
    all_w1 = np.asarray([r["macro_f1_w1"] for r in paired], dtype=np.float64)
    all_r = np.asarray([r["mean_random_f1"] for r in paired], dtype=np.float64)
    overlaps = []
    for k in range(n_folds):
        a, b = panel_rows[2 * k:2 * k + 2]
        shared = a["genes"] & b["genes"]
        overlaps.append({"fold": k, "shared_genes": len(shared), "w0_only_genes": k_genes - len(shared),
                         "w1_only_genes": k_genes - len(shared), "jaccard": len(shared) / len(a["genes"] | b["genes"]),
                         "w0_genes": a["genes"], "w1_genes": b["genes"]})
    stats = dict(n_individuals=len(paired), n_folds=n_folds, n_random_draws_per_fold=n_draws,
                 paired_delta=float(all_d.mean()), median_delta=float(np.median(all_d)),
                 mean_macro_f1_w0=float(all_w0.mean()), mean_macro_f1_w1=float(all_w1.mean()),
                 random_reference_mean=float(all_r.mean()), w0_minus_random=float(all_w0.mean() - all_r.mean()),
                 w1_minus_random=float(all_w1.mean() - all_r.mean()),
                 per_individual_positive=int(sign["positive"]), per_individual_negative=int(sign["negative"]),
                 per_individual_zero=int(sign["zero"]),
                 w0_exceeds_all_random_fold_means=all(q["macro_f1_w0"] > q["random_max_draw"] for q in byfold),
                 individuals_w0_exceeds_every_random_draw=int(beats_all),
                 panels=[{k2: (sorted(v) if isinstance(v, set) else v) for k2, v in r.items()} for r in panel_rows],
                 source_sha256=sources)
    return paired, byfold, overlaps, stats


def csv_text(rows, columns) -> str:
    buf = io.StringIO(newline="")
    w = csv.DictWriter(buf, fieldnames=columns, extrasaction="ignore", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
    return buf.getvalue()


def tables(paired, byfold, overlaps) -> dict[str, bytes]:
    clean = [{k: v for k, v in q.items() if k != "draws"} for q in byfold]
    return {
        "table_01_individuals_anonymized.csv": csv_text(paired, list(paired[0].keys())).encode(),
        "table_02_fold_summary.csv": csv_text(clean, list(clean[0].keys())).encode(),
        "table_03_panel_overlap.csv": csv_text([{k: v for k, v in r.items() if k not in ("w0_genes", "w1_genes")} for r in overlaps],
                                               ["fold", "shared_genes", "w0_only_genes", "w1_only_genes", "jaccard"]).encode(),
        "table_04_selected_gene_membership.csv": csv_text(
            [{"fold": q["fold"], "gene": g, "in_w0": int(g in q["w0_genes"]), "in_w1": int(g in q["w1_genes"])}
             for q in overlaps for g in sorted(q["w0_genes"] | q["w1_genes"])], ["fold", "gene", "in_w0", "in_w1"]).encode(),
    }


def figures(paired, byfold, overlaps, stats, k_genes: int) -> dict[str, bytes]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    n_folds = len(byfold)
    out = {}
    per_fold = [sum(1 for r in paired if r["fold"] == k) for k in range(n_folds)]
    starts = np.cumsum([0] + per_fold)

    fig, ax = plt.subplots(figsize=(9, 5.1))
    for k in range(n_folds):
        vals = sorted(r["paired_delta"] for r in paired if r["fold"] == k)
        xx = np.arange(starts[k] + 1, starts[k] + len(vals) + 1)
        ax.scatter(xx, vals, s=20, alpha=0.8, color=FOLD_COLORS[k % 4], label=f"Fold {k} (n={len(vals)})")
        if k != n_folds - 1:
            ax.axvline(starts[k + 1] + 0.5, lw=0.7, ls=":", color=GREY, alpha=0.8)
    mean = stats["paired_delta"]
    ax.axhline(0, lw=1, ls="--", color="#444444", label="No difference")
    ax.axhline(mean, lw=1.7, color="#111111", label=f"Donor-equal mean {mean:+.4f}")
    ax.set(xlim=(0, starts[-1] + 1), xlabel="Anonymous donor rank within each fold", ylabel="Macro-F1(w=0) − Macro-F1(w=1)",
           title=f"Paired donor-level differences (exploratory; {len(paired)} donors)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=3, fontsize=8.5, frameon=False)
    fig.tight_layout()
    out[FIGURE_NAMES[0]] = fig

    fig, ax = plt.subplots(figsize=(8.7, 4.9))
    xx = np.arange(n_folds)
    for q in byfold:
        ax.scatter([xx[q["fold"]]] * len(q["draws"]), q["draws"], marker="_", s=110, color=GREY, alpha=0.75, linewidths=1.3,
                   label=f"Random draw (fold mean, {len(q['draws'])} per fold)" if q["fold"] == 0 else None)
    ax.scatter(xx, [q["mean_random_f1"] for q in byfold], marker="D", s=46, color="#555b63", zorder=3, label="Random-reference mean")
    ax.scatter(xx - 0.09, [q["macro_f1_w0"] for q in byfold], marker="o", s=58, color=W0C, zorder=4, label="w=0 (cell class)")
    ax.scatter(xx + 0.09, [q["macro_f1_w1"] for q in byfold], marker="^", s=64, color=W1C, zorder=4, label="w=1 (cell class + donor)")
    ax.set(xlabel="Donor-held-out CV fold", ylabel="Donor-equal mean macro-F1", xlim=(-0.5, n_folds - 0.5), xticks=xx,
           xticklabels=[f"Fold {k}" for k in range(n_folds)], title="Compact panels versus expression-matched random reference")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2, fontsize=8.5, frameon=False)
    fig.tight_layout()
    out[FIGURE_NAMES[1]] = fig

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.bar([f"Fold {k}" for k in range(n_folds)], [o["jaccard"] for o in overlaps], color=W0C)
    for k, q in enumerate(overlaps):
        ax.text(k, q["jaccard"] + 0.015, f"{q['shared_genes']}/{k_genes} shared\nJ={q['jaccard']:.3f}", ha="center", va="bottom",
                fontsize=9)
    ax.axhline(0.8, ls="--", lw=1.2, color="#444444")
    ax.text(n_folds - 0.55, 0.81, "Proposal's informative-overlap threshold (J = 0.8)", ha="right", va="bottom", fontsize=8.5,
            color="#333333")
    ax.set(ylim=(0, 0.92), ylabel="Jaccard(w=0, w=1)", title=f"Gene-panel overlap by held-out fold ({k_genes} genes each)")
    fig.tight_layout()
    out[FIGURE_NAMES[2]] = fig

    blobs = {}
    for name, f in out.items():
        buf = io.BytesIO()
        f.savefig(buf, format="png", dpi=180, bbox_inches="tight")
        plt.close(f)
        blobs[name] = buf.getvalue()
    return blobs


EXPLORATORY_LABEL = ("Exploratory numerical sensitivity checks (G3 fixed-alpha tolerance refits, G4 identity-permutation "
                     "sanity check, H1 fallback diagnosis). They are not confirmatory gates and do not replace "
                     "pre-specified controls such as a decoy task, a displacement control or fresh-seed repetitions.")


def parameters(cfg) -> dict:
    return {"classes": cfg.classes, "n_folds": cfg.n_folds, "k": cfg.k, "n_hvg": int(cfg["hvg"]["n_top_genes"]),
            "min_cells_per_class": int(cfg["cohort"]["min_cells_per_class"]),
            "min_evaluable_classes": int(cfg["cohort"]["min_evaluable_classes"]),
            "withholding_fraction": str(cfg["withholding"]["fraction"]), "selector_cap": int(cfg["subsample"]["selector_cap"]),
            "classifier_cap": int(cfg["subsample"]["classifier_cap"]), "random_draws": int(cfg["random_reference"]["draws"]),
            "random_bins": int(cfg["random_reference"]["bins"]), "stability_tol": float(cfg["stability"]["tol"]),
            "g3_threshold": float(cfg["stability"]["g3_threshold"]), "g4_threshold": float(cfg["stability"]["g4_threshold"])}


def dataset_record(run_dir: Path) -> dict:
    """Dataset dimensions and cohort counts from the most recent passing validation record of the run."""
    vd = run_dir / "validation"
    reps = sorted(vd.glob("validation_*.json")) if vd.exists() else []
    for p in reversed(reps):
        r = read_json(p)
        if r["status"] != "PASS":
            continue
        c = {x["id"]: x for x in r["checks"]}
        shape = c.get("C3", {}).get("x", {}).get("shape")
        obs = c.get("C4", {}).get("observed", {})
        return {"n_cells_total": shape[0] if shape else None, "n_genes_total": shape[1] if shape else None,
                "n_eligible_individuals": obs.get("n_eligible"), "n_analysis_cells": obs.get("n_analysis_cells"),
                "eligible_by_site": obs.get("eligible_by_site"), "source": f"validation/{p.name}"}
    return {}


def stability_summary(run_dir: Path) -> dict:
    sd = run_dir / "stability"
    out = {}
    if not sd.exists():
        return out
    for p in sorted(sd.glob("g3_fold*_*.json")):
        r = read_json(p)
        out[p.stem] = {"status": r["status"], "fixed_alpha": r["fixed_alpha"], "support_size": r["support_size"],
                       "reference_exact64": r["reference_exact64"],
                       "support_jaccard": r["support_jaccard"], "top64_jaccard": r["top64_jaccard_if_applicable"],
                       "n_iter": r["fit"]["n_iter"], "converged": r["fit"]["converged"]}
    for name in ("g4_state", "h1_state"):
        p = sd / f"{name}.json"
        if p.exists():
            r = read_json(p)
            res = r.get("result", r.get("conclusion", {}))
            entry = {"status": r["status"], "trace": [(x["alpha"], x["nnz"], x["n_iter"]) for x in r["records"]],
                     "result": {k: v for k, v in res.items() if k not in ("selected", "reference")}}
            if "selected" in res:
                a, b = set(res["selected"]), set(res["reference"])
                entry["result"].update(shared_with_reference=len(a & b), union_with_reference=len(a | b))
            out[name.replace("_state", "")] = entry
    return out


# ---------------------------------------------------------------------- regenerated report
_WORDS = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine", 10: "ten"}


def _w(n: int, cap: bool = False) -> str:
    w = _WORDS.get(int(n), str(n))
    return w.capitalize() if cap else w


def _pow10(x: float) -> str:
    import math
    e = math.log10(x)
    sup = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
    return f"10{str(int(round(e))).translate(sup)}" if abs(e - round(e)) < 1e-12 else f"{x:g}"


def _s(x: float, nd: int) -> str:
    """Signed number with a typographic minus."""
    return f"{x:+.{nd}f}".replace("-", "−")


def check_template_has_no_numbers(template: str) -> None:
    """The report template may contain no result-like numbers: every quantity must be a placeholder.
    Allowed: version labels, the reference list, section/figure/table numbers and small integers in prose."""
    import re
    body = template.split("## References")[0]
    body = re.sub(r"\{\{\w+\}\}", "", body)
    body = re.sub(r"v\d+\.\d+", "", body)
    body = re.sub(r"(?m)SHA-256|\(\d{4}\)|^#+ \d+(\.\d+)*", "", body)
    bad = re.findall(r"\d+[.,]\d+|\d{3,}", body)
    require(not bad, f"report template contains hard-coded numbers: {bad[:10]}")


def render_report(template: str, results: dict, comparison: dict | None, results_sha: str, build_date: str) -> str:
    import re
    s, folds, ov = results["summary"], results["fold_summary"], results["panel_overlap"]
    st, par, ds = results["stability"]["checks"], results["parameters"], results.get("dataset", {})
    check_template_has_no_numbers(template)
    K = par["k"]
    rows = "\n".join(f"| {q['fold']} | {q['donors']} | {q['macro_f1_w0']:.4f} | {q['macro_f1_w1']:.4f} | "
                     f"{q['mean_random_f1']:.4f} | {_s(q['paired_delta'], 4)} |" for q in folds)
    fd = [q["paired_delta"] for q in folds]
    jac = [o["jaccard"] for o in ov]
    shared = [o["shared_genes"] for o in ov]
    g3 = {k: v for k, v in st.items() if k.startswith("g3_")}
    stable = [k for k, v in g3.items() if v["status"] in ("STABLE_AT_FIXED_ALPHA", "FALLBACK_TOP64_STABLE_AT_FIXED_ALPHA")]
    fb = [p for p in s["panels"] if not p["exact64"]]
    if g3:
        extra = f", including the truncated top-{K} set of the fallback panel" if any(
            g3[k]["status"].startswith("FALLBACK") for k in stable) else ""
        g3_sentence = (f"Fixed-alpha, tighter-tolerance (tol {_pow10(par['stability_tol'])}) cold refits recovered the same {K} selected genes for "
                       f"{len(stable)} of {len(g3)} panels{extra}.")
        det = "; ".join(f"fold {k[7]} w={k[-1]}: {v['support_size']} active genes, support Jaccard {v['support_jaccard']:.4f}"
                        + (f", top-{K} Jaccard {v['top64_jaccard']:.4f}" if not v.get("reference_exact64", True) else "")
                        for k, v in sorted(g3.items()))
        g3_detail = f"{g3_sentence} Per panel — {det}."
    else:
        g3_sentence = g3_detail = "The fixed-alpha tolerance checks were not run."
    g4 = st.get("g4")
    g4_sentence = (f"The identity-permutation sanity check (fold 0) selected a panel with Jaccard "
                   f"{g4['result']['shared_with_reference']}/{g4['result']['union_with_reference']} = {g4['result']['jaccard']:.4f} "
                   f"relative to the w=0 panel ({g4['result']['diagnostic_verdict']} at the {g4['result']['threshold']:.2f} "
                   f"threshold; {len(g4['trace'])} cold fits).") if g4 and g4.get("result") else "The permutation check was not run."
    h1 = st.get("h1")
    h1_sentence = (f"At {len(h1['trace'])} fixed alphas around the fallback alpha, cold fits had "
                   f"{', '.join(str(t[1]) for t in h1['trace'])} active genes; "
                   f"{'an exact-' + str(K) + ' solution occurred' if h1['result'].get('exact_64_found') else 'no exact-' + str(K) + ' solution occurred'} at these alphas."
                   ) if h1 else ""
    fallback_sentence = " ".join(f"The fold {p['fold']} w={p['arm'][-1]} panel used the fallback: {K + p['discarded']} active genes "
                                 f"reduced to the {K} largest coefficient-row norms." for p in fb) or f"All panels had exactly {K} active genes."
    if comparison:
        c = comparison["counts"]
        n = sum(c.values())
        repro_sentence = (f"all {n} compared components agreed bit-identically with the historical records."
                          if c.get("EXACT", 0) == n else
                          f"of {n} compared components, {c.get('EXACT', 0)} agreed bit-identically, "
                          f"{c.get('NUMERICAL_AGREEMENT', 0)} within tolerance and {c.get('DISCREPANCY', 0)} showed discrepancies.")
        repro_detail = comparison.get("narrative", repro_sentence)
    else:
        repro_sentence = repro_detail = "no historical comparison was available for this build."
    values = {
        "build_date": build_date, "results_sha": results_sha, "n_individuals": s["n_individuals"],
        "n_cells": f"{sum(q['test_cells'] for q in folds):,}", "delta_4": _s(s["paired_delta"], 4), "delta_7": _s(s["paired_delta"], 7),
        "median_4": _s(s["median_delta"], 4), "pos": s["per_individual_positive"], "neg": s["per_individual_negative"],
        "zero": s["per_individual_zero"], "fold_delta_min": _s(min(fd), 4), "fold_delta_max": _s(max(fd), 4),
        "k": K, "n_folds": _w(len(folds)), "n_folds_cap": _w(len(folds), True), "n_panels": _w(2 * len(folds)),
        "n_panels_cap": _w(2 * len(folds), True), "n_classes": _w(len(par["classes"])),
        "min_classes_word": _w(par["min_evaluable_classes"]),
        "class_list": ", ".join(par["classes"]), "min_cells": par["min_cells_per_class"], "min_classes": par["min_evaluable_classes"],
        "n_hvg": f"{par['n_hvg']:,}", "per_fold_individuals": "–".join(sorted({str(q['donors']) for q in folds})),
        "dataset_cells": f"{ds['n_cells_total']:,}" if ds.get("n_cells_total") else "(not recorded)",
        "dataset_genes": f"{ds['n_genes_total']:,}" if ds.get("n_genes_total") else "(not recorded)",
        "exploratory_label": results["stability"]["label"],
        "jac_min": f"{min(jac):.2f}", "jac_max": f"{max(jac):.2f}", "jac_mean": f"{sum(jac) / len(jac):.2f}",
        "shared_min": min(shared), "shared_max": max(shared), "changed_min": K - max(shared), "changed_max": K - min(shared),
        "n_draws": s["n_random_draws_per_fold"], "n_random_models": s["n_random_draws_per_fold"] * len(folds),
        "n_classifiers": _w(2 * len(folds)),
        "w0_minus_random_3": _s(s["w0_minus_random"], 3), "w1_minus_random_3": _s(s["w1_minus_random"], 3),
        "w0_4": f"{s['mean_macro_f1_w0']:.4f}", "w1_4": f"{s['mean_macro_f1_w1']:.4f}", "random_4": f"{s['random_reference_mean']:.4f}",
        "w0_6": f"{s['mean_macro_f1_w0']:.6f}", "w1_6": f"{s['mean_macro_f1_w1']:.6f}", "random_6": f"{s['random_reference_mean']:.6f}",
        "beats_fold_means": "yes, in all folds" if s["w0_exceeds_all_random_fold_means"] else "not in every fold",
        "beats_all": s["individuals_w0_exceeds_every_random_draw"], "table_rows": rows, "g3_sentence": g3_sentence,
        "g3_detail": g3_detail, "g4_sentence": g4_sentence, "h1_sentence": h1_sentence, "fallback_sentence": fallback_sentence,
        "repro_sentence": repro_sentence, "repro_detail": repro_detail,
    }
    out = re.sub(r"\{\{(\w+)\}\}", lambda m: str(values[m.group(1)]), template)
    require("{{" not in out, "unresolved report placeholder")
    return out


def build_results(run_dir: Path, out_dir: Path, cfg, *, comparison: dict | None = None, overwrite: bool = False,
                  extra_files: dict[str, bytes] | None = None) -> dict:
    paired, byfold, overlaps, stats = collect(run_dir, cfg)
    files = tables(paired, byfold, overlaps)
    files.update(figures(paired, byfold, overlaps, stats, cfg.k))
    results = {"schema": "compact-gene-panels/results/v1", "analysis": "descriptive, exploratory, one repetition",
               "config_sha256": cfg.sha256, "summary": stats,
               "fold_summary": [{k: v for k, v in q.items()} for q in byfold],
               "panel_overlap": [{k: (sorted(v) if isinstance(v, set) else v) for k, v in r.items()} for r in overlaps],
               "parameters": parameters(cfg), "dataset": dataset_record(run_dir),
               "stability": {"label": EXPLORATORY_LABEL, "checks": stability_summary(run_dir)}}
    if comparison is not None:
        results["comparison_with_historical_reference"] = {"overall": comparison["overall"], "counts": comparison["counts"],
                                                           "tolerances": comparison["tolerances"]}
    files["results.json"] = json_bytes(results)
    files.update(extra_files or {})
    for name, blob in sorted(files.items()):
        write_bytes(out_dir / name, blob, overwrite=overwrite)
    manifest = "".join(f"{sha256_bytes(blob)}  {name}\n" for name, blob in sorted(files.items()))
    write_bytes(out_dir / "SHA256SUMS.txt", manifest.encode(), overwrite=overwrite)
    write_bytes(out_dir / "build_record.json", json_bytes({"built_utc": utc_now(), "environment": environment(),
                                                           "files": {n: sha256_bytes(b) for n, b in sorted(files.items())}}),
                overwrite=True)
    return {"files": sorted(files), "summary": stats}


def render_readme_block(results: dict, comparison: dict | None) -> str:
    s, ov = results["summary"], results["panel_overlap"]
    jac = [o["jaccard"] for o in ov]
    rows = [("Held-out individuals (all folds)", f"{s['n_individuals']}"),
            ("Mean macro-F1, w=0 (cell class)", f"{s['mean_macro_f1_w0']:.6f}"),
            ("Mean macro-F1, w=1 (cell class + identity)", f"{s['mean_macro_f1_w1']:.6f}"),
            ("Mean macro-F1, expression-matched random panels", f"{s['random_reference_mean']:.6f}"),
            ("Donor-equal mean paired Δ = F1(w=0) − F1(w=1)", _s(s["paired_delta"], 7)),
            ("Individuals with Δ > 0 / < 0 / = 0", f"{s['per_individual_positive']} / {s['per_individual_negative']} / {s['per_individual_zero']}"),
            ("Within-fold panel Jaccard(w=0, w=1)", f"{min(jac):.3f}–{max(jac):.3f}")]
    out = ["*Generated by `cgp rebuild-results` from `results/results.json`; descriptive, one cross-validation repetition.*", "",
           "| Quantity | Value |", "|---|---|"] + [f"| {a} | {b} |" for a, b in rows]
    out += ["", "A negative Δ means w=1 scored slightly higher; it does not indicate benefit, equivalence or non-inferiority. "
            "Identity supervision changed many selected genes with little change in held-out performance.", "",
            "![Paired differences](figures/figure_01_paired_donor_deltas.png)"]
    if comparison:
        c = comparison["counts"]
        out += ["", f"**Reproduction (Level A) against the historical records:** overall `{comparison['overall']}` — "
                f"{c.get('EXACT', 0)} components bit-identical, {c.get('NUMERICAL_AGREEMENT', 0)} within tolerance, "
                f"{c.get('DISCREPANCY', 0)} discrepancies ({comparison['narrative']})"]
    return "\n".join(out)
