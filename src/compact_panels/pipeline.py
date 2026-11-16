"""Stage orchestration with isolated, checkpointed, write-once outputs.

Run directory layout (``<run_dir>``)::

    run_manifest.json                 configuration hash, environment, input hashes
    logs/<stage>.log
    fold<k>/plan.json                 list sizes and SHA-256 of every cell / individual list
    fold<k>/preprocessing.npz         training-only HVG indices, genes, means, SDs
    fold<k>/checkpoint_<arm>.npz      alpha-search state after every fit (resumable)
    fold<k>/panel_<arm>.json          selected 64-gene panel + alpha trace summary
    fold<k>/model_<arm>.npz/.json     logistic-regression parameters + diagnostics
    fold<k>/training_complete.json
    fold<k>/random/manifest.json, model_<dd>.npz, training_complete.json
    fold<k>/evaluation_attempt.json   written BEFORE any held-out expression is read
    fold<k>/evaluation.json           per-individual scores (main arms + random draws)
    stability/*.json

Training stages read expression of fitting cells only; scoring reads held-out cells
once and refuses to run twice.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from . import classifier as clf
from . import data, evaluation, identity, partition, preprocessing, random_reference, selection, stability, targets
from .io_utils import IntegrityError, read_json, require, sha256_file, sha256_lines, write_json, write_npz
from .provenance import environment, maxrss, utc_now

ARMS = ("w0", "w1")


def input_spec(cfg) -> dict:
    """Optional ``inputs.toml`` next to the study configuration: provenance of derivable inputs."""
    import tomllib
    p = cfg.path.parent / "inputs.toml"
    if not p.exists():
        return {}
    with open(p, "rb") as f:
        return tomllib.load(f)


class Logger:
    def __init__(self, path: Path | None):
        self.path = path
        if path:
            path.parent.mkdir(parents=True, exist_ok=True)

    def __call__(self, msg: str) -> None:
        line = f"{utc_now()} {msg}"
        print(line, flush=True)
        if self.path:
            with open(self.path, "a", encoding="utf-8") as f:
                f.write(line + "\n")


class Study:
    """Lazy access to metadata, cohort, fold assignment and fold plans for one run."""

    def __init__(self, cfg, paths, log=None, verify_inputs: bool = True):
        self.cfg, self.paths = cfg, paths
        self.log = log or Logger(None)
        self._meta = self._cohort = self._assign = None
        self._plans: dict[int, partition.FoldPlan] = {}
        self.input_hashes: dict[str, str] = {}
        if verify_inputs:
            self.verify_inputs()

    # ------------------------------------------------------------------ inputs
    def verify_inputs(self) -> dict:
        """File SHA-256 must equal the pin; a re-derived identity mapping (``cgp derive-inputs``) is
        accepted instead when its (patient_id, individual_id) pair hash equals the pin in inputs.toml."""
        exp = self.cfg.inputs
        derived = input_spec(self.cfg)
        for key, p in (("h5ad", self.paths.h5ad), ("identity_mapping", self.paths.identity_mapping),
                       ("qc_exclusion", self.paths.qc_exclusion)):
            require(p.is_file(), f"input file missing: {key} ({p.name})")
            h = sha256_file(p)
            self.input_hashes[key] = h
            want = exp.get(f"{key}_sha256")
            if want and h != want and key == "identity_mapping" and derived.get("identity_mapping", {}).get("pairs_sha256"):
                pairs = identity.mapping_pairs_sha(identity.load_identity_mapping(p))
                require(pairs == derived["identity_mapping"]["pairs_sha256"],
                        f"identity mapping differs from the verified pairs (pair SHA-256 {pairs})")
                self.input_hashes["identity_mapping_pairs"] = pairs
                continue
            if want:
                require(h == want, f"SHA-256 mismatch for {key}: expected {want}, found {h}")
        return self.input_hashes

    def metadata(self):
        if self._meta is None:
            self._meta = data.load_metadata(self.paths.h5ad)
        return self._meta

    def cohort(self) -> identity.Cohort:
        if self._cohort is None:
            obs, _, _ = self.metadata()
            mapping = identity.load_identity_mapping(self.paths.identity_mapping)
            excl = identity.load_cell_list(self.paths.qc_exclusion)
            self._cohort = identity.build_cohort(obs, mapping, excl, self.cfg)
        return self._cohort

    def assignment(self) -> partition.FoldAssignment:
        if self._assign is None:
            self._assign = partition.assign_folds(self.cohort(), self.cfg)
        return self._assign

    def plan(self, fold: int) -> partition.FoldPlan:
        if fold not in self._plans:
            self._plans[fold] = partition.build_fold_plan(self.cohort(), self.assignment(), fold, self.cfg)
        return self._plans[fold]

    # ------------------------------------------------------------------ run directory
    def fold_dir(self, fold: int) -> Path:
        return self.paths.run_dir / f"fold{fold}"

    def ensure_run_manifest(self) -> dict:
        path = self.paths.run_dir / "run_manifest.json"
        if path.exists():
            m = read_json(path)
            require(m["config_sha256"] == self.cfg.sha256,
                    "run directory was created with a different study configuration; use a new run_dir")
            require(m["input_sha256"] == self.input_hashes, "run directory was created from different inputs")
            return m
        m = {"created_utc": utc_now(), "config_sha256": self.cfg.sha256, "config": self.cfg.raw,
             "config_file": self.cfg.path.name, "input_sha256": self.input_hashes, "environment": environment(),
             "paths": self.paths.as_public_dict()}
        write_json(path, m)
        return m

    def write_plan(self, fold: int) -> dict:
        plan = self.plan(fold)
        rec = {"fold": fold, "lists": plan.provenance(), "sha": plan.provenance_sha(), "seeds": plan.seeds,
               "fold_assignment_sha256": self.assignment().canonical_sha256,
               "convention": self.cfg.convention(fold)}
        path = self.fold_dir(fold) / "plan.json"
        if path.exists():
            require(read_json(path)["sha"] == rec["sha"], f"fold {fold}: plan changed since the run started")
        else:
            write_json(path, rec)
        return rec

    # ------------------------------------------------------------------ expression
    def read_rows(self, names: list[str], allowed: set[str]):
        obs, _, _ = self.metadata()
        return data.read_expression(self.paths.h5ad, names, obs.index, allowed=allowed)

    def preprocess(self, fold: int):
        """Training-only HVG + standardization; writes preprocessing.npz once, then verifies it."""
        plan = self.plan(fold)
        conv = self.cfg.convention(fold)
        blk = self.read_rows(plan.fitting, allowed=set(plan.fitting))
        hv = self.cfg["hvg"]
        pre, F = preprocessing.fit_preprocessing(blk.matrix, blk.var_names, int(hv["n_top_genes"]), hv["flavor"],
                                                 conv["hvg_convention"])
        del blk
        path = self.fold_dir(fold) / "preprocessing.npz"
        if path.exists():
            with np.load(path, allow_pickle=False) as z:
                require(np.array_equal(pre.candidate_var_idx, z["candidate_var_idx"])
                        and np.array_equal(pre.genes, z["genes"])
                        and np.array_equal(pre.means, z["means"]) and np.array_equal(pre.sds, z["sds"]),
                        f"fold {fold}: preprocessing differs from the saved checkpoint")
        else:
            write_npz(path, candidate_var_idx=pre.candidate_var_idx, genes=pre.genes, means=pre.means, sds=pre.sds)
        return pre, F

    # ------------------------------------------------------------------ stages
    def run_panels(self, fold: int, max_new_fits: int) -> bool:
        require(max_new_fits > 0, "max_new_fits must be positive")
        dest = self.fold_dir(fold)
        require(not (dest / "evaluation_attempt.json").exists(), f"fold {fold} already opened for evaluation")
        self.ensure_run_manifest()
        prov = self.write_plan(fold)
        if all((dest / f"panel_{a}.json").exists() for a in ARMS):
            self.log(f"fold {fold}: both panels already complete")
            return True
        plan = self.plan(fold)
        self.log(f"fold {fold}: {len(plan.fitting)} fitting, {len(plan.withheld)} withheld, "
                 f"{len(plan.selector)} selector, {len(plan.classifier)} classifier, {len(plan.test)} held-out cells")
        pre, F = self.preprocess(fold)
        pos = preprocessing.positions(plan.fitting, plan.selector)
        X = preprocessing.standardize_rows(F, pos, pre.means, pre.sds, fortran=True)
        del F
        obs = self.cohort().obs.loc[plan.selector]
        ctype = self.cfg["data"]["cell_type_column"]
        provenance = {"plan_sha": prov["sha"], "preprocessing_sha": sha256_file(dest / "preprocessing.npz"),
                      "config_sha256": self.cfg.sha256}
        remaining = max_new_fits
        for arm in ARMS:
            if not (dest / f"panel_{arm}.json").exists() and remaining <= 0:
                self.log(f"fold {fold}: fit budget used; rerun the same command to resume")
                return False
            Y = targets.targets_for_arm(arm, obs[ctype].to_numpy(), obs.individual_id.to_numpy(), plan.train_people,
                                        self.cfg.classes)
            used = selection.run_search(X, Y, arm=arm, fold=fold, genes=pre.genes, cfg=self.cfg, dest=dest,
                                        max_new_fits=max(remaining, 1), provenance=provenance, log=self.log)
            remaining -= used
            if not (dest / f"panel_{arm}.json").exists():
                self.log(f"fold {fold}: fit budget used; rerun the same command to resume")
                return False
        return True

    def _classifier_matrix(self, fold: int):
        plan = self.plan(fold)
        dest = self.fold_dir(fold)
        pre_path = dest / "preprocessing.npz"
        require(pre_path.exists(), f"fold {fold}: preprocessing missing")
        pre, F = self.preprocess(fold)  # verifies the checkpoint
        pos = preprocessing.positions(plan.fitting, plan.classifier)
        X = preprocessing.standardize_rows(F, pos, pre.means, pre.sds, fortran=False)
        del F
        y = self.cohort().obs.loc[plan.classifier][self.cfg["data"]["cell_type_column"]].astype(str).to_numpy()
        require(set(y) == set(self.cfg.classes), "classifier subsample lacks a class")
        return pre, X, y

    def run_classifiers(self, fold: int) -> None:
        dest = self.fold_dir(fold)
        require(not (dest / "evaluation_attempt.json").exists(), f"fold {fold} already opened for evaluation")
        if (dest / "training_complete.json").exists():
            self.log(f"fold {fold}: classifiers already trained")
            return
        for arm in ARMS:
            require((dest / f"panel_{arm}.json").exists(), f"fold {fold}: panel {arm} missing; run panels first")
        prov = self.write_plan(fold)
        pre, X, y = self._classifier_matrix(fold)
        for arm in ARMS:
            panel = selection.load_panel(dest, arm)
            require(panel["all_converged"] and panel["provenance"]["preprocessing_sha"] == sha256_file(dest / "preprocessing.npz"),
                    f"fold {fold}/{arm}: panel provenance mismatch")
            ii = np.asarray(panel["panel_idx"], dtype=np.int64)
            res = clf.train_logistic(X[:, ii], y, self.cfg)
            require(res["converged"] and set(res["classes"]) == set(self.cfg.classes), f"fold {fold}/{arm}: classifier did not converge")
            write_npz(dest / f"model_{arm}.npz", coef=res["coef"], intercept=res["intercept"], classes=res["classes"],
                      panel_idx=ii, gene_names=pre.genes[ii].astype(str))
            write_json(dest / f"model_{arm}.json", {
                "fold": fold, "arm": arm, "converged": True, "n_iter": res["n_iter"], "seconds": res["seconds"],
                "convergence_warnings": res["convergence_warnings"], "other_warnings": res["other_warnings"],
                "hyperparameters": res["hyperparameters"], "n_cells": int(len(y)), "n_genes": int(len(ii)),
                "class_counts": {c: int(np.sum(y == c)) for c in self.cfg.classes},
                "model_sha256": sha256_file(dest / f"model_{arm}.npz"), "panel_sha256": selection.panel_sha(dest, arm),
                "plan_sha": prov["sha"], "peak_rss_bytes": maxrss()})
            self.log(f"fold {fold}/{arm}: classifier trained ({res['n_iter']} iterations, {res['seconds']:.2f}s)")
        write_json(dest / "training_complete.json", {
            "status": "TRAINED_NO_HELD_OUT_EXPRESSION_READ", "fold": fold, "plan_sha": prov["sha"],
            "preprocessing_sha": sha256_file(dest / "preprocessing.npz"),
            "panels": {a: selection.panel_sha(dest, a) for a in ARMS},
            "models": {a: sha256_file(dest / f"model_{a}.npz") for a in ARMS}, "completed_utc": utc_now()})

    def run_random(self, fold: int) -> None:
        dest = self.fold_dir(fold) / "random"
        fdir = self.fold_dir(fold)
        require(not (fdir / "evaluation_attempt.json").exists(), f"fold {fold} already opened for evaluation")
        if (dest / "training_complete.json").exists():
            self.log(f"fold {fold}: random-reference models already trained")
            return
        require((fdir / "panel_w0.json").exists(), f"fold {fold}: w0 panel missing")
        prov = self.write_plan(fold)
        rr = self.cfg["random_reference"]
        pre, X, y = self._classifier_matrix(fold)
        w0 = selection.load_panel(fdir, "w0")
        panels, counts = random_reference.random_panels(pre.means, w0["panel_idx"], fold, n_genes=len(pre.genes),
                                                        k=self.cfg.k, bins=int(rr["bins"]), draws=int(rr["draws"]),
                                                        seed=int(rr["seed"]), salt=int(rr["salt"]))
        manifest = {"fold": fold, "draws": int(rr["draws"]), "k": self.cfg.k, "bins": int(rr["bins"]),
                    "seed_rule": [int(rr["seed"]), fold, "draw", int(rr["salt"])], "panels": panels, "bin_counts": counts,
                    "preprocessing_sha": sha256_file(fdir / "preprocessing.npz"),
                    "w0_panel_sha": selection.panel_sha(fdir, "w0"), "plan_sha": prov["sha"]}
        mpath = dest / "manifest.json"
        if mpath.exists():
            require(read_json(mpath) == json.loads(json.dumps(manifest)), f"fold {fold}: random manifest mismatch")
        else:
            write_json(mpath, manifest)
        for d, panel in enumerate(panels):
            path = dest / f"model_{d:02d}.npz"
            if path.exists():
                with np.load(path, allow_pickle=False) as z:
                    require(np.array_equal(z["cols"], panel), f"fold {fold}: random draw {d} gene mismatch")
                continue
            res = clf.train_logistic(X[:, np.asarray(panel)], y, self.cfg)
            require(res["converged"] and res["coef"].shape == (len(self.cfg.classes), self.cfg.k),
                    f"fold {fold}: random draw {d} did not converge")
            write_npz(path, coef=res["coef"], intercept=res["intercept"], classes=res["classes"],
                      cols=np.asarray(panel, dtype=np.int64), n_iter=np.array(res["n_iter"]))
        write_json(dest / "training_complete.json", {
            "fold": fold, "manifest_sha": sha256_file(mpath),
            "models_sha": {f"model_{d:02d}": sha256_file(dest / f"model_{d:02d}.npz") for d in range(len(panels))},
            "held_out_expression_read": False, "completed_utc": utc_now()})
        self.log(f"fold {fold}: {len(panels)} random-reference classifiers trained")

    def run_score(self, fold: int) -> dict:
        """Read held-out expression ONCE and score w0, w1 and every random draw."""
        fdir, rdir = self.fold_dir(fold), self.fold_dir(fold) / "random"
        attempt, out = fdir / "evaluation_attempt.json", fdir / "evaluation.json"
        require(not attempt.exists() and not out.exists(), f"fold {fold}: held-out set already opened; refusing to rescore")
        tc = fdir / "training_complete.json"
        require(tc.exists() and (rdir / "training_complete.json").exists(), f"fold {fold}: training incomplete")
        rec = read_json(tc)
        prov = self.write_plan(fold)
        require(rec["plan_sha"] == prov["sha"], "fold membership changed since training")
        pre_path = fdir / "preprocessing.npz"
        require(sha256_file(pre_path) == rec["preprocessing_sha"], "preprocessing changed since training")
        with np.load(pre_path, allow_pickle=False) as z:
            cand, genes, mu, sd = z["candidate_var_idx"].copy(), z["genes"].astype(str), z["means"].copy(), z["sds"].copy()
        models = {}
        for arm in ARMS:
            mp = fdir / f"model_{arm}.npz"
            require(sha256_file(mp) == rec["models"][arm], f"fold {fold}/{arm}: model changed since training")
            with np.load(mp, allow_pickle=False) as z:
                models[arm] = (z["coef"].copy(), z["intercept"].copy(), z["classes"].astype(str), z["panel_idx"].copy())
        rrec = read_json(rdir / "training_complete.json")
        rmodels = []
        for key in sorted(rrec["models_sha"]):
            mp = rdir / f"{key}.npz"
            require(sha256_file(mp) == rrec["models_sha"][key], f"fold {fold}: random model {key} changed")
            with np.load(mp, allow_pickle=False) as z:
                rmodels.append((z["coef"].copy(), z["intercept"].copy(), z["classes"].astype(str), z["cols"].copy()))
        plan = self.plan(fold)
        write_json(attempt, {"fold": fold, "started_utc": utc_now(), "test_cells_sha": sha256_lines(plan.test),
                             "training_complete_sha": sha256_file(tc), "config_sha256": self.cfg.sha256})
        needed = sorted({int(i) for m in models.values() for i in m[3]} | {int(i) for m in rmodels for i in m[3]})
        index = {g: j for j, g in enumerate(needed)}
        blk = self.read_rows(plan.test, allowed=set(plan.test))
        require(np.array_equal(np.asarray(blk.var_names[cand], dtype=str), genes), "held-out gene columns differ")
        Xall = blk.matrix[:, cand[needed]].toarray().astype(np.float64)
        del blk
        Xall = (Xall - mu[needed]) / sd[needed]
        require(np.all(np.isfinite(Xall)), "non-finite held-out values")
        meta = self.cohort().obs.loc[plan.test]
        ctype = self.cfg["data"]["cell_type_column"]
        truth = meta[ctype].astype(str).to_numpy()
        donors = meta.individual_id.astype(str).to_numpy()
        ev = self.cfg["evaluation"]
        mc, mk = int(ev["min_cells_per_class"]), int(ev["min_evaluable_classes"])
        preds = {arm: clf.predict(c, b, cl, Xall[:, [index[int(i)] for i in ii]]) for arm, (c, b, cl, ii) in models.items()}
        rows = evaluation.score_donors(truth, donors, preds, plan.test_people, self.cfg.classes, mc, mk)
        for r in rows:
            r["paired_delta"] = r["macro_f1_w0"] - r["macro_f1_w1"]
        draws = []
        for c, b, cl, cols in rmodels:
            p = clf.predict(c, b, cl, Xall[:, [index[int(i)] for i in cols]])
            draws.append([r["macro_f1_r"] for r in evaluation.score_donors(truth, donors, {"r": p}, plan.test_people,
                                                                            self.cfg.classes, mc, mk)])
        V = np.asarray(draws, dtype=np.float64)  # draws x donors
        for i, r in enumerate(rows):
            r["f1_by_random_draw"] = V[:, i].tolist()
            r["mean_random_f1"] = float(V[:, i].mean())
        require(sum(r["n_test_cells"] for r in rows) == len(plan.test), "scoring coverage error")
        result = {"fold": fold, "n_individuals": len(rows), "n_test_cells": len(plan.test),
                  "mean_macro_f1_w0": float(np.mean([r["macro_f1_w0"] for r in rows])),
                  "mean_macro_f1_w1": float(np.mean([r["macro_f1_w1"] for r in rows])),
                  "paired_delta": float(np.mean([r["paired_delta"] for r in rows])),
                  "mean_random_f1": float(V.mean()), "random_mean_by_draw": V.mean(axis=1).tolist(),
                  "per_individual": rows, "training_complete_sha": sha256_file(tc),
                  "random_training_complete_sha": sha256_file(rdir / "training_complete.json"),
                  "completed_utc": utc_now()}
        write_json(out, result)
        self.log(f"fold {fold}: scored once; {len(rows)} individuals; delta={result['paired_delta']:+.6f}")
        return result

    # ------------------------------------------------------------------ stability
    def selector_matrices(self, fold: int):
        plan = self.plan(fold)
        pre, F = self.preprocess(fold)
        pos = preprocessing.positions(plan.fitting, plan.selector)
        X = preprocessing.standardize_rows(F, pos, pre.means, pre.sds, fortran=True)
        del F
        obs = self.cohort().obs.loc[plan.selector]
        cls = obs[self.cfg["data"]["cell_type_column"]].to_numpy()
        return pre, X, cls, obs.individual_id.to_numpy(), plan.train_people

    def run_g3(self, fold: int, arm: str) -> dict:
        st = self.cfg["stability"]
        out = self.paths.run_dir / "stability" / f"g3_fold{fold}_{arm}.json"
        require(not out.exists(), f"{out.name} exists (write-once)")
        fdir = self.fold_dir(fold)
        panel = selection.load_panel(fdir, arm)
        pre, X, cls, ind, people = self.selector_matrices(fold)
        Y = targets.targets_for_arm(arm, cls, ind, people, self.cfg.classes)
        alpha = float(panel["alpha_final"])
        norms, fit = stability.cold_fit(X, Y, alpha, tol=float(st["tol"]), max_iter=int(st["max_iter"]))
        support = np.flatnonzero(norms > 0).tolist()
        k = self.cfg.k
        top = stability.top_k_positive(norms, k) if len(support) >= k else None
        status, raw_j, top_j = stability.judge_fixed_alpha(panel["panel_idx"], panel["exact64"], support, top,
                                                           float(st["g3_threshold"]))
        if not fit["converged"]:
            status = "NONCONVERGED_REVIEW_REQUIRED"
        res = {"check": "G3_fixed_alpha", "fold": fold, "arm": arm, "status": status, "fixed_alpha": alpha,
               "reference_exact64": bool(panel["exact64"]), "fit": fit, "support_size": len(support),
               "support_jaccard": raw_j, "top64_jaccard_if_applicable": top_j,
               "tight_support_genes": pre.genes[support].tolist(),
               "tight_top64_genes": pre.genes[top].tolist() if top is not None else None,
               "panel_sha256": selection.panel_sha(fdir, arm)}
        write_json(out, res)
        self.log(f"G3 fold{fold}/{arm}: {status} support={len(support)} J={raw_j:.4f} top={top_j}")
        return res

    def run_h1(self) -> dict:
        st = self.cfg["stability"]
        fold, arm = int(st["h1_fold"]), str(st["h1_arm"])
        out = self.paths.run_dir / "stability" / "h1_state.json"
        panel = selection.load_panel(self.fold_dir(fold), arm)
        pre, X, cls, ind, people = self.selector_matrices(fold)
        Y = targets.targets_for_arm(arm, cls, ind, people, self.cfg.classes)
        state = read_json(out) if out.exists() else {"check": "H1_fixed_alphas", "fold": fold, "arm": arm, "records": [],
                                                       "status": "RUNNING"}
        done = {float(r["alpha"]) for r in state["records"]}
        k = self.cfg.k
        for alpha in st["h1_alphas"]:
            if float(alpha) in done:
                continue
            norms, rec = stability.cold_fit(X, Y, float(alpha), tol=float(st["tol"]), max_iter=int(st["max_iter"]))
            if not rec["converged"]:
                state["records"].append(rec); state["status"] = "STOP_REVIEW"
                write_json(out, state, overwrite=True)
                raise IntegrityError("H1 cold fit did not converge")
            idx = stability.top_k_positive(norms, k) if rec["nnz"] >= k else []
            rec.update(selected_gene_names=pre.genes[np.flatnonzero(norms > 0)].tolist(),
                       top64_gene_names=pre.genes[idx].tolist() if idx else [],
                       jaccard_panel_if_64=stability.jaccard(idx, panel["panel_idx"]) if idx else None)
            state["records"].append(rec)
            write_json(out, state, overwrite=True)
            self.log(f"H1 alpha={alpha} active={rec['nnz']} iter={rec['n_iter']}")
        nnz = [r["nnz"] for r in state["records"]]
        state.update(status="COMPLETE_DIAGNOSTIC", conclusion={"n_active_across_tested_alphas": nnz,
                                                                "uniform_support": len(set(nnz)) == 1,
                                                                "exact_64_found": k in nnz})
        write_json(out, state, overwrite=True)
        return state

    def run_g4(self, max_new_fits: int = 50) -> dict:
        st = self.cfg["stability"]
        sel = self.cfg["selector"]
        fold = int(st["g4_fold"])
        out = self.paths.run_dir / "stability" / "g4_state.json"
        w0 = selection.load_panel(self.fold_dir(fold), "w0")
        pre, X, cls, ind, people = self.selector_matrices(fold)
        perm = np.random.default_rng(int(st["g4_seed"])).permutation(len(ind))
        permuted = ind[perm]
        require(sorted(permuted.tolist()) == sorted(ind.tolist()), "permutation changed donor counts")
        Y = targets.make_y(cls, permuted, people, self.cfg.classes, 1.0)
        am = float(selection.alpha_max(X, Y))
        k = self.cfg.k
        if out.exists():
            state = read_json(out)
            require(np.isclose(state["alpha_max"], am, rtol=1e-12), "G4 state provenance changed")
            if state["status"] in ("COMPLETE", "STOP_REVIEW"):
                return state
        else:
            state = {"check": "G4_identity_permutation", "fold": fold, "status": "RUNNING", "alpha_max": am,
                     "low": am * float(sel["lower_frac"]), "high": am, "next_alpha": am * float(sel["first_alpha_frac"]),
                     "records": [], "best_alpha_over_k": None, "best_norms_over_k": None,
                     "permutation_seed": int(st["g4_seed"])}
        limit = int(sel["max_steps"]) + int(sel["extra_steps"])
        while len(state["records"]) < limit and max_new_fits > 0:
            alpha = float(state["next_alpha"])
            norms, rec = stability.cold_fit(X, Y, alpha, tol=float(st["tol"]), max_iter=int(st["max_iter"]))
            rec["step"] = len(state["records"]) + 1
            state["records"].append(rec)
            max_new_fits -= 1
            if not rec["converged"]:
                state["status"] = "STOP_REVIEW"
                write_json(out, state, overwrite=True)
                raise IntegrityError("G4 fit did not converge")
            state, idx, exact = stability.advance_cold_search(
                state, alpha=alpha, rec=rec, norms=norms, k=k, max_steps=int(sel["max_steps"]),
                extra_steps=int(sel["extra_steps"]), trigger=int(sel["fallback_trigger"]), lower_frac=float(sel["lower_frac"]))
            self.log(f"G4 step {rec['step']}: alpha={alpha:.8g} active={rec['nnz']} iter={rec['n_iter']}")
            if idx is not None:
                jac = stability.jaccard(idx, w0["panel_idx"])
                state["status"] = "COMPLETE"
                state["result"] = {"selected": pre.genes[idx].tolist(), "reference": w0["panel_genes"], "jaccard": jac,
                                   "threshold": float(st["g4_threshold"]),
                                   "diagnostic_verdict": "PASS" if jac >= float(st["g4_threshold"]) else "FAIL",
                                   "exact_64": exact, "selected_alpha": alpha if exact else state["best_alpha_over_k"],
                                   "n_fits": rec["step"]}
            write_json(out, state, overwrite=True)
            if state["status"] != "RUNNING":
                break
        return state
