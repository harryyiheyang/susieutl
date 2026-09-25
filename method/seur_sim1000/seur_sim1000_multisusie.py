import argparse
import csv
import os
import pickle
import sys
import time

import numpy as np


ap = argparse.ArgumentParser()
ap.add_argument("--rep-root", required=True)
ap.add_argument("--env-root", required=True)
args = ap.parse_args()
root = os.path.abspath(args.rep_root)
env = os.path.abspath(args.env_root)
ad = os.path.join(root, "adapter")
od = os.path.join(root, "methods", "multisusie")
if not os.path.isdir(ad):
    raise RuntimeError("adapter output is missing")
if os.path.exists(od):
    raise RuntimeError("refusing to overwrite MultiSuSiE output: " + od)
os.makedirs(od)
src = os.path.join(env, "vendor", "MultiSuSiE", "src")
if not os.path.isdir(src):
    raise RuntimeError("official MultiSuSiE source is missing")
sys.path.insert(0, src)
from MultiSuSiE.susiepy_ss import multisusie_rss, susie_get_pip

with open(os.path.join(ad, "adapter_manifest.tsv"), newline="", encoding="utf-8") as f:
    manifest = list(csv.DictReader(f, delimiter="\t"))
if not manifest or len({x["input_md5"] for x in manifest}) != 1:
    raise RuntimeError("adapter input identity is invalid")
with open(os.path.join(ad, "variants.tsv"), newline="", encoding="utf-8") as f:
    variants = list(csv.DictReader(f, delimiter="\t"))
vn = [x["SNP"] for x in variants]
if vn != [f"v{i:04d}" for i in range(1, 501)]:
    raise RuntimeError("adapter variant order is not canonical")
an = ["EUR", "EAS", "AFR"]
n = [1000000, 150000, 150000]
R = []
z = []
for a in an:
    rv = np.fromfile(os.path.join(ad, "multi", f"R_{a}.f64le"), dtype="<f8")
    zv = np.fromfile(os.path.join(ad, "multi", f"z_{a}.f64le"), dtype="<f8")
    if rv.size != 500 * 500 or zv.size != 500:
        raise RuntimeError("adapter binary dimensions are invalid for " + a)
    Ra = rv.reshape((500, 500), order="F")
    if not np.all(np.isfinite(Ra)) or not np.all(np.isfinite(zv)):
        raise RuntimeError("adapter has non-finite values for " + a)
    if np.max(np.abs(Ra - Ra.T)) > 1e-10 or np.max(np.abs(np.diag(Ra) - 1)) > 1e-10:
        raise RuntimeError("adapter LD semantic check failed for " + a)
    R.append(Ra)
    z.append(zv)

t0 = time.perf_counter()
fit = multisusie_rss(R_list=R, population_sizes=n, z_list=z, L=5)
elapsed = time.perf_counter() - t0
pip = np.asarray(susie_get_pip(fit), dtype=np.float64)
if pip.shape != (500,) or not np.all(np.isfinite(pip)) or np.any(pip < 0) or np.any(pip > 1):
    raise RuntimeError("MultiSuSiE native PIP is invalid")
ranks = np.array([1 + np.sum(pip > x) for x in pip], dtype=int)
with open(os.path.join(od, "pip.tsv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["variant", "variant_row", "pip", "rank"])
    for j in range(500):
        w.writerow([vn[j], j + 1, format(pip[j], ".17g"), int(ranks[j])])

coef = np.asarray(fit.coef, dtype=np.float64)
coef_sd = np.asarray(fit.coef_sd, dtype=np.float64)
moment_ok = (coef.shape == (3, 500) and coef_sd.shape == (3, 500) and
             np.all(np.isfinite(coef)) and np.all(np.isfinite(coef_sd)) and
             np.all(coef_sd >= 0))
metric_status = "available_native_posterior_moments" if moment_ok else "unavailable_native_moment_contract"
if moment_ok:
    from math import erfc, sqrt
    with open(os.path.join(od, "effects.tsv"), "w", newline="", encoding="utf-8") as ef, \
         open(os.path.join(od, "lfsr.tsv"), "w", newline="", encoding="utf-8") as lf:
        ew = csv.writer(ef, delimiter="\t", lineterminator="\n")
        lw = csv.writer(lf, delimiter="\t", lineterminator="\n")
        ew.writerow(["variant", "variant_row", "ancestry", "causal", "mean", "variance", "metric_status"])
        lw.writerow(["variant", "variant_row", "ancestry", "causal", "value", "call_0.10", "call_0.05", "metric_status"])
        for t, a in enumerate(an):
            for j in range(500):
                sd = coef_sd[t, j]
                value = 1.0 if sd == 0 and coef[t, j] == 0 else (0.0 if sd == 0 else 0.5 * erfc(abs(coef[t, j] / sd) / sqrt(2.0)))
                causal = int(j + 1 in (209, 280, 386))
                ew.writerow([vn[j], j + 1, a, causal, format(coef[t, j], ".17g"), format(sd * sd, ".17g"), metric_status])
                lw.writerow([vn[j], j + 1, a, causal, format(value, ".17g"), int(value <= .10), int(value <= .05), metric_status])
else:
    for name in ("effects.tsv", "lfsr.tsv"):
        with open(os.path.join(od, name), "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(["metric_status"])
            w.writerow([metric_status])

with open(os.path.join(od, "cs.tsv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["cs_id", "variants", "size"])
    sets = fit.sets.get("cs", {}) if isinstance(fit.sets, dict) else {}
    for key, value in sets.items():
        idx = np.asarray(value, dtype=int).reshape(-1)
        w.writerow([key, ",".join(vn[j] for j in idx), int(idx.size)])
with open(os.path.join(od, "status.tsv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["method", "status", "elapsed_sec", "input_md5", "version", "commit", "V_policy", "metric_status", "call"])
    w.writerow(["multisusie", "success", format(elapsed, ".17g"), manifest[0]["input_md5"], "source",
                "10351ebbbc0ed4c891c455000eb6ec7800d0b717", "native_auto_default", metric_status,
                "multisusie_rss(R_list,population_sizes,z_list,L=5)"])
with open(os.path.join(od, "fit.pkl"), "wb") as f:
    pickle.dump(fit, f, protocol=pickle.HIGHEST_PROTOCOL)
