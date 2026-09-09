"""Section 1b precompute: taper-count (K) sweep on AF.SKRH-XV.BAND at fixed NW.

Classical Mspec at K in {3,7,13,21,29,45} (same bandwidth product NW, so a fair
apples-to-apples comparison at fixed resolution), each in a fresh subprocess so
runtime and peak memory are clean isolated numbers. Shows: (i) how fast the
coherence estimate's roughness stabilizes with K, and (ii) that runtime and
memory keep climbing well past the point roughness stops improving.

FastMspec's own fused single point (auto K, from driver_sec1's output) is
overlaid by the notebook, not recomputed here.

Output: nb3_sec1b_ksweep.json -- [{k_taps, roughness, runtime_s, peak_mem_mb, coh_num}]
"""
import json, os, subprocess, sys
import numpy as np

D = "/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch"
sys.path.insert(0, D + "/python")
from dispcurve_pick.gvib_loader import build_pair_matched_data

GVIB = "/scratch/tolugboj_lab/Prj5_HarnomicRFTraces/para_prepross/ADAMA_gvib.h5"
STALIST = D + "/ADAMA_stalist.csv"
OUT = D + "/report_figures"
PY = sys.executable
WORKER = os.path.join(os.path.dirname(os.path.abspath(__file__)), "worker.py")

NW_PROD = 0.001 * 10801
K_VALUES = [3, 7, 13, 21, 29, 45]


def roughness(est, f, lo=0.0, hi=0.4):
    band = (f >= lo) & (f < hi)
    e = est[band]
    d2 = e[2:] - 2 * e[1:-1] + e[:-2]
    return float(np.mean(np.abs(d2)) / np.sqrt(np.mean(e ** 2)))


print("Building AF.SKRH-XV.BAND from gvib ...", flush=True)
pmd = build_pair_matched_data(GVIB, "AF.SKRH", "XV.BAND", STALIST)
s1_npy, s2_npy = OUT + "/_skrhb_s1.npy", OUT + "/_skrhb_s2.npy"
np.save(s1_npy, pmd.S1_data_mat)
np.save(s2_npy, pmd.S2_data_mat)

rows = []
for k in K_VALUES:
    out_npz = f"{OUT}/_kw_{k}.npz"
    print(f"  Mspec K={k} ...", flush=True)
    p = subprocess.run([PY, WORKER, s1_npy, s2_npy, "Mspec", str(NW_PROD), out_npz, str(k)],
                       capture_output=True, text=True)
    if p.returncode:
        sys.stdout.write(p.stderr[-2000:])
        raise SystemExit(f"K={k} failed")
    m = json.loads([l for l in p.stdout.strip().splitlines() if l.startswith("{")][-1])
    d = np.load(out_npz)
    m["roughness"] = round(roughness(d["coh"], d["freq"]), 4)
    m["k_taps"] = k
    rows.append(m)
    print(f"     rough={m['roughness']}  {m['runtime_s']}s  {m['peak_mem_mb']} MB", flush=True)
    os.remove(out_npz)

with open(OUT + "/nb3_sec1b_ksweep.json", "w") as f:
    json.dump(rows, f, indent=2)
for npy in (s1_npy, s2_npy):
    os.remove(npy)
print("DONE ->", OUT + "/nb3_sec1b_ksweep.json", flush=True)
print(json.dumps(rows, indent=2))
