"""Command-line interface: ``cgp <command>``.

All commands take ``--config`` (study configuration, default ``configs/study.toml``) and
``--paths`` (paths file, default ``$CGP_PATHS`` or ``configs/paths.local.toml``). Commands
are idempotent: rerunning a completed stage verifies and skips it; an interrupted panel
search resumes from its last checkpoint.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import tomllib
from pathlib import Path

from . import __version__
from .config import load_paths, load_study
from .io_utils import IntegrityError, read_json, write_json
from .pipeline import ARMS, Logger, Study


def _default(path: str) -> str | None:
    return path if Path(path).exists() else None


def _common(p: argparse.ArgumentParser) -> None:
    p.add_argument("--config", default="configs/study.toml", help="study configuration TOML")
    p.add_argument("--paths", default=None, help="paths TOML (default: $CGP_PATHS or configs/paths.local.toml)")
    p.add_argument("--run-dir", default=None, help="override the run directory")
    p.add_argument("--refs", default="configs/historical_references.toml", help="historical reference pins")
    p.add_argument("--acceptance", default="configs/acceptance.toml", help="comparison tolerances")


def _setup(a, stage: str, verify_inputs: bool = True):
    cfg = load_study(a.config)
    paths = load_paths(a.paths or os.environ.get("CGP_PATHS") or _default("configs/paths.local.toml"), run_dir=a.run_dir)
    log = Logger(paths.run_dir / "logs" / f"{stage}.log")
    return cfg, paths, Study(cfg, paths, log, verify_inputs=verify_inputs), log


def _tol(a) -> dict:
    with open(a.acceptance, "rb") as f:
        return tomllib.load(f)["tolerances"]


def _folds(a, cfg) -> list[int]:
    return list(range(cfg.n_folds)) if a.fold is None else [a.fold]


def cmd_validate(a) -> int:
    from .validate import run_validation
    cfg = load_study(a.config)
    paths = load_paths(a.paths or os.environ.get("CGP_PATHS") or _default("configs/paths.local.toml"), run_dir=a.run_dir)
    refs = Path(a.refs) if a.refs and Path(a.refs).exists() else None
    tol = _tol(a) if Path(a.acceptance).exists() else None
    rep = run_validation(cfg, paths, refs_path=refs, acceptance=tol, deep=a.deep)
    print(f"VALIDATION {rep['status']}: {sum(c['passed'] for c in rep['checks'])}/{len(rep['checks'])} checks passed; "
          f"models fitted: 0; held-out expression read: no")
    return 0 if rep["status"] == "PASS" else 1


def cmd_plan(a) -> int:
    cfg, paths, st, log = _setup(a, "plan")
    st.ensure_run_manifest()
    for k in _folds(a, cfg):
        rec = st.write_plan(k)
        log(f"fold {k}: " + ", ".join(f"{n}={v['size']}" for n, v in rec["lists"].items()))
    return 0


def cmd_panels(a) -> int:
    cfg, paths, st, log = _setup(a, "panels")
    ok = True
    for k in _folds(a, cfg):
        ok &= st.run_panels(k, a.max_new_fits)
    return 0 if ok else 3


def cmd_classifiers(a) -> int:
    cfg, paths, st, log = _setup(a, "classifiers")
    for k in _folds(a, cfg):
        st.run_classifiers(k)
    return 0


def cmd_random(a) -> int:
    cfg, paths, st, log = _setup(a, "random")
    for k in _folds(a, cfg):
        st.run_random(k)
    return 0


def cmd_score(a) -> int:
    if not a.confirm_held_out_access:
        print("Scoring reads held-out expression exactly once per fold. Re-run with --confirm-held-out-access.")
        return 2
    cfg, paths, st, log = _setup(a, "score")
    for k in _folds(a, cfg):
        if (st.fold_dir(k) / "evaluation.json").exists():
            log(f"fold {k}: already scored; not repeated")
            continue
        st.run_score(k)
    return 0


def cmd_stability(a) -> int:
    cfg, paths, st, log = _setup(a, "stability")
    if a.check == "g3":
        folds = _folds(a, cfg)
        arms = ARMS if a.arm is None else (a.arm,)
        for k in folds:
            for arm in arms:
                if (paths.run_dir / "stability" / f"g3_fold{k}_{arm}.json").exists():
                    log(f"G3 fold{k}/{arm}: already done")
                    continue
                st.run_g3(k, arm)
    elif a.check == "g4":
        st.run_g4(a.max_new_fits)
    elif a.check == "h1":
        if cfg["stability"].get("h1_alphas"):
            st.run_h1()
        else:
            log("H1 not configured for this study configuration")
    return 0


def cmd_compare(a) -> int:
    from .comparison import build_comparison
    from .reference import Legacy
    cfg, paths, st, log = _setup(a, "compare", verify_inputs=False)
    if not paths.legacy_root:
        print("compare requires legacy_root in the paths file")
        return 2
    res = build_comparison(paths.run_dir, Legacy(paths.legacy_root, Path(a.refs)), cfg.n_folds, _tol(a))
    print(json.dumps({"overall": res["overall"], "counts": res["counts"]}, indent=1))
    return 0 if res["overall"] != "DISCREPANCY" else 4


def cmd_report(a) -> int:
    from .reporting import build_results
    cfg, paths, st, log = _setup(a, "report", verify_inputs=False)
    comp = paths.run_dir / "comparisons" / "latest.json"
    out = Path(a.out) if a.out else paths.run_dir / "results"
    res = build_results(paths.run_dir, out, cfg, comparison=read_json(comp) if comp.exists() else None, overwrite=a.overwrite)
    s = res["summary"]
    log(f"results written to {out.name}/: delta={s['paired_delta']:+.7f} w0={s['mean_macro_f1_w0']:.6f} "
        f"w1={s['mean_macro_f1_w1']:.6f} random={s['random_reference_mean']:.6f}")
    return 0


def cmd_rebuild_results(a) -> int:
    """Rebuild results/, figures/ and the regenerated report from completed run outputs."""
    import shutil
    from .comparison import build_comparison
    from .io_utils import json_bytes, sha256_file, write_bytes
    from .provenance import utc_now
    from .reference import Legacy
    from .reporting import FIGURE_NAMES, build_results, render_report
    cfg, paths, st, log = _setup(a, "rebuild", verify_inputs=False)
    repo = Path(a.repo_root).resolve()
    comparison, extra = None, {}
    if paths.legacy_root and a.refs and Path(a.refs).exists():
        comparison = build_comparison(paths.run_dir, Legacy(paths.legacy_root, Path(a.refs)), cfg.n_folds, _tol(a))
        extra["comparison_with_historical_reference.json"] = json_bytes(comparison)
    res_dir = repo / "results"
    if comparison is not None:
        import tempfile
        from .comparison import compare_results_package
        from .reference import worst
        with tempfile.TemporaryDirectory() as tmp:          # first pass: package to compare with the historical delivery
            build_results(paths.run_dir, Path(tmp), cfg)
            pkg = compare_results_package(Path(tmp), Legacy(paths.legacy_root, Path(a.refs)), _tol(a))
        comparison["results_package"] = pkg
        sts = [v["status"] for v in pkg.values()]
        comparison["by_component"]["results package (tables, summary values, figures)"] = {x: sts.count(x) for x in set(sts)}
        for x in sts:
            comparison["counts"][x] = comparison["counts"].get(x, 0) + 1
        comparison["overall"] = worst([comparison["overall"]] + sts)
        comparison["narrative"] = "; ".join(f"{n}: {', '.join(f'{k} {s}' for s, k in sorted(c.items()))}"
                                            for n, c in comparison["by_component"].items()) + "."
        extra["comparison_with_historical_reference.json"] = json_bytes(comparison)
    build_results(paths.run_dir, res_dir, cfg, comparison=comparison, overwrite=True, extra_files=extra)
    (repo / "figures").mkdir(exist_ok=True)
    for name in FIGURE_NAMES:
        shutil.copyfile(res_dir / name, repo / "figures" / name)
    results = read_json(res_dir / "results.json")
    text = render_report((repo / "report" / "report_template.md").read_text(encoding="utf-8"), results, comparison,
                         sha256_file(res_dir / "results.json"), utc_now()[:10])
    from .reporting import render_readme_block
    readme = repo / "README.md"
    if readme.exists():
        import re as _re
        txt = readme.read_text(encoding="utf-8")
        block = render_readme_block(results, comparison)
        new = _re.sub(r"<!-- RESULTS:START -->.*?<!-- RESULTS:END -->",
                      lambda m: f"<!-- RESULTS:START -->\n{block}\n<!-- RESULTS:END -->", txt, flags=_re.S)
        write_bytes(readme, new.encode("utf-8"), overwrite=True)
    report_path = repo / "report" / f"compact_gene_panels_report_v{a.report_version}.md"
    write_bytes(report_path, text.encode("utf-8"), overwrite=True)
    files = sorted([p for d in ("results", "figures") for p in (repo / d).iterdir() if p.is_file()] + [report_path])
    manifest = "".join(f"{sha256_file(p)}  {p.relative_to(repo)}\n" for p in files)
    write_bytes(repo / "RELEASE_SHA256SUMS.txt", manifest.encode(), overwrite=True)
    log(f"rebuilt {len(files)} release files; comparison overall: {comparison['overall'] if comparison else 'not run'}")
    return 0 if not comparison or comparison["overall"] != "DISCREPANCY" else 4


def cmd_reproduce(a) -> int:
    """Sequential driver: plan -> panels -> classifiers -> random [-> score -> stability -> report]."""
    cfg, paths, st, log = _setup(a, "reproduce")
    st.ensure_run_manifest()
    folds = _folds(a, cfg)
    for k in folds:
        while not st.run_panels(k, a.max_new_fits):
            pass
        st.run_classifiers(k)
        st.run_random(k)
    if a.through in ("score", "all"):
        if not a.confirm_held_out_access:
            print("Stopping before scoring: pass --confirm-held-out-access to score held-out individuals.")
            return 0
        for k in folds:
            if not (st.fold_dir(k) / "evaluation.json").exists():
                st.run_score(k)
    if a.through == "all":
        for k in folds:
            for arm in ARMS:
                if not (paths.run_dir / "stability" / f"g3_fold{k}_{arm}.json").exists():
                    st.run_g3(k, arm)
        if int(cfg["stability"]["g4_fold"]) in folds:
            st.run_g4()
        if cfg["stability"].get("h1_alphas") and int(cfg["stability"]["h1_fold"]) in folds:
            st.run_h1()
    return 0


def cmd_status(a) -> int:
    cfg = load_study(a.config)
    paths = load_paths(a.paths or os.environ.get("CGP_PATHS") or _default("configs/paths.local.toml"), run_dir=a.run_dir)
    for k in range(cfg.n_folds):
        d = paths.run_dir / f"fold{k}"
        items = ["plan.json", "preprocessing.npz", "panel_w0.json", "panel_w1.json", "training_complete.json",
                 "random/training_complete.json", "evaluation.json"]
        print(f"fold {k}: " + "  ".join(f"{i.split('/')[0].split('.')[0]}={'Y' if (d / i).exists() else '-'}" for i in items))
    return 0


def cmd_derive_inputs(a) -> int:
    """Re-derive the identity mapping (public SDRF) and the QC exclusion list (metadata rule)."""
    from . import data, identity
    from .io_utils import sha256_file, sha256_lines, write_bytes, write_lines
    from .pipeline import input_spec
    cfg = load_study(a.config)
    spec = input_spec(cfg)
    h5ad, out = Path(a.h5ad), Path(a.out)
    want = spec.get("h5ad", {}).get("sha256")
    got = sha256_file(h5ad)
    if want and got != want:
        print(f"H5AD SHA-256 mismatch: {got}")
        return 5
    obs, _, _ = data.load_metadata(h5ad)
    mapping = identity.mapping_from_sdrf(a.sdrf, obs[cfg["data"]["patient_column"]].astype(str))
    pairs = identity.mapping_pairs_sha(mapping)
    qc = identity.qc_exclusion_from_rule(obs, spec["qc_exclusion"]["rule"])
    ok_map = pairs == spec["identity_mapping"]["pairs_sha256"]
    ok_qc = sha256_lines(qc) == spec["qc_exclusion"]["sha256"]
    write_bytes(out / "identity_mapping.tsv", mapping.to_csv(sep="\t", index=False, lineterminator="\n").encode())
    write_lines(out / "qc_exclusion_cells.txt", qc)
    print(f"identity mapping: {len(mapping)} rows, pair SHA-256 {'MATCHES' if ok_map else 'DIFFERS FROM'} the pin")
    print(f"QC exclusion list: {len(qc)} cells, SHA-256 {'MATCHES' if ok_qc else 'DIFFERS FROM'} the pin")
    return 0 if ok_map and ok_qc else 5


def cmd_synthetic_demo(a) -> int:
    """End-to-end demonstration on generated synthetic data (no real data needed)."""
    from .synthetic import make_synthetic
    out = Path(a.out).resolve()
    if out.exists() and any(out.iterdir()):
        print(f"{out} is not empty; choose a new directory")
        return 2
    info = make_synthetic(out / "data")
    (out / "paths.toml").write_text(
        "[paths]\nh5ad = \"data/synthetic.h5ad\"\nidentity_mapping = \"data/identity_mapping.tsv\"\n"
        "qc_exclusion = \"data/qc_exclusion_cells.txt\"\nrun_dir = \"run\"\n")
    ns = argparse.Namespace(config=a.config, paths=str(out / "paths.toml"), run_dir=None, refs=None, acceptance=a.acceptance,
                            fold=None, max_new_fits=50, through="all", confirm_held_out_access=True, deep=False)
    from .validate import run_validation
    cfg = load_study(a.config)
    paths = load_paths(ns.paths)
    rep = run_validation(cfg, paths, refs_path=None, acceptance=None, deep=True)
    if rep["status"] != "PASS":
        print("synthetic validation failed")
        return 1
    cmd_reproduce(ns)
    ns.out = str(out / "results")
    ns.overwrite = False
    cmd_report(ns)
    print(f"SYNTHETIC DEMO COMPLETE: {info['n_cells']} cells x {info['n_genes']} genes; results in {out / 'results'}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="cgp", description=__doc__.splitlines()[0])
    ap.add_argument("--version", action="version", version=f"compact-gene-panels {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)
    specs = {
        "validate": (cmd_validate, "read-only input validation and dry run (no fitting)"),
        "plan": (cmd_plan, "write fold plans (metadata only)"),
        "panels": (cmd_panels, "fit w0/w1 MultiTaskLasso panels (resumable)"),
        "classifiers": (cmd_classifiers, "train the panel classifiers"),
        "random": (cmd_random, "train the expression-matched random-reference classifiers"),
        "score": (cmd_score, "score held-out individuals once per fold"),
        "stability": (cmd_stability, "training-only numerical stability checks"),
        "compare": (cmd_compare, "compare the run with the historical reference artifacts"),
        "report": (cmd_report, "build tables, figures, results.json and manifest from the run"),
        "reproduce": (cmd_reproduce, "sequential end-to-end driver"),
        "rebuild-results": (cmd_rebuild_results, "rebuild results/, figures/ and the report from completed run outputs"),
        "status": (cmd_status, "show stage completion per fold"),
        "derive-inputs": (cmd_derive_inputs, "re-derive the identity mapping and QC exclusion list from public sources"),
        "synthetic-demo": (cmd_synthetic_demo, "end-to-end demonstration on synthetic data"),
    }
    for name, (fn, help_) in specs.items():
        p = sub.add_parser(name, help=help_)
        _common(p)
        p.set_defaults(fn=fn)
        if name in ("plan", "panels", "classifiers", "random", "score", "stability", "reproduce"):
            p.add_argument("--fold", type=int, default=None)
        if name in ("panels", "stability", "reproduce"):
            p.add_argument("--max-new-fits", type=int, default=50)
        if name in ("score", "reproduce"):
            p.add_argument("--confirm-held-out-access", action="store_true")
        if name == "validate":
            p.add_argument("--deep", action="store_true", help="also replay training-only preprocessing (reads fitting rows)")
        if name == "stability":
            p.add_argument("check", choices=("g3", "g4", "h1"))
            p.add_argument("--arm", choices=ARMS, default=None)
        if name == "reproduce":
            p.add_argument("--through", choices=("train", "score", "all"), default="train")
        if name in ("report", "synthetic-demo"):
            p.add_argument("--out", default=None)
        if name == "report":
            p.add_argument("--overwrite", action="store_true")
        if name == "rebuild-results":
            p.add_argument("--repo-root", default=".")
            p.add_argument("--report-version", default="2.0")
        if name == "derive-inputs":
            p.add_argument("--h5ad", required=True)
            p.add_argument("--sdrf", required=True, help="E-MTAB-10026 SDRF file downloaded from ArrayExpress")
            p.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    if a.cmd == "synthetic-demo" and a.config == "configs/study.toml":
        a.config = "configs/synthetic.toml"
    from .selection import SearchStop
    try:
        return int(a.fn(a) or 0)
    except (IntegrityError, SearchStop) as e:
        print(str(e), file=sys.stderr)
        return 5


if __name__ == "__main__":
    sys.exit(main())
