"""Section 1 + 1c precompute for Notebook 3.

Section 1: AF.SKRH-XV.BAND, built from ADAMA_gvib.h5, all four techniques at
their *production* configs -- single-taper; MspecBestK and FastMspec at
Wband=0.001 (NW~=10.8); classical Mspec at the batch pipeline's fixed
NW=100 / K=80 (work_unit.py's NW_MSPEC / K_MSPEC). Coherence spectrum +
isolated runtime + peak memory for each, one fresh process per run.

Section 1c: FastMspec vs MspecBestK swept over NW in {10.8, 18, 24, 28, 40}.
10.8-28 stays at or below SKRH-BAND's real resolution ceiling
(NW_high = N*c_min/(4R) ~= 25-30 for a physical Love-wave c_min ~3 km/s);
40 is deliberately past it -- flagged as a cost probe, not a usable estimate.
Shows FastMspec's cost (driven by the bounded transition-region count r)
staying flat while MspecBestK's (driven by K, which grows with NW) climbs.

Every run is wrapped so one failure records an error row and the job continues.

Outputs (small, committed to the repo):
  nb3_sec1_skrh.npz        -- freq + 4 coherence curves (Section 1, production configs)
  nb3_sec1_nwsweep.npz     -- freq + FastMspec/MspecBestK coherence at each swept NW
  nb3_sec1_metrics.json    -- [{section, technique, nw, nw_label, config, runtime_s, peak_mem_mb, taper_size, coh_num, error}]
"""
import json, os, subprocess, sys, traceback
import numpy as np

D = "/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch"
sys.path.insert(0, D + "/python")
from dispcurve_pick.gvib_loader import build_pair_matched_data

GVIB = "/scratch/tolugboj_lab/Prj5_HarnomicRFTraces/para_prepross/ADAMA_gvib.h5"
STALIST = D + "/ADAMA_stalist.csv"
OUT = D + "/report_figures"
PY = sys.executable
WORKER = os.path.join(os.path.dirname(os.path.abspath(__file__)), "worker.py")
os.makedirs(OUT, exist_ok=True)

N = 10801
NW_PROD = 0.001 * N            # 10.801  (Wband=0.001)
NW_SWEEP = [18.0, 24.0, 28.0, 40.0]   # 10.8 reused from Section 1; 28 ~= ceiling; 40 past it
SAFE = {"single-taper": "single_taper", "MspecBestK": "MspecBestK", "FastMspec": "FastMspec", "Mspec": "Mspec"}

print(f"Building AF.SKRH-XV.BAND matched data from {GVIB} ...", flush=True)
pmd = build_pair_matched_data(GVIB, "AF.SKRH", "XV.BAND", STALIST)
print(f"  dist={pmd.dist_km:.1f} km, shape={pmd.S1_data_mat.shape}, {len(pmd.days_used)} days", flush=True)
s1_npy, s2_npy = OUT + "/_skrh_s1.npy", OUT + "/_skrh_s2.npy"
np.save(s1_npy, pmd.S1_data_mat)
np.save(s2_npy, pmd.S2_data_mat)

metrics = []
curves = {}       # (label, technique) -> coh array
freq_holder = {}  # "freq" -> array (kept out of `curves` so its keys stay uniform 2-tuples)


def run(technique, nw, tag, section, nw_label, config, k_taps=None):
    out_npz = f"{OUT}/_w_{technique}_{tag}.npz"
    argv = [PY, WORKER, s1_npy, s2_npy, technique, str(nw), out_npz]
    if k_taps is not None:
        argv.append(str(k_taps))
    print(f"  [{section}] {technique} {nw_label} ({config}) ...", flush=True)
    row = dict(section=section, technique=technique, nw=round(nw, 3), nw_label=nw_label,
               config=config, runtime_s=None, peak_mem_mb=None, taper_size=None, coh_num=None, error=None)
    try:
        p = subprocess.run(argv, capture_output=True, text=True, timeout=4 * 3600)
        if p.returncode != 0:
            row["error"] = (p.stderr or "")[-600:]
            print(f"     FAILED: {row['error'][-200:]}", flush=True)
        else:
            m = json.loads([l for l in p.stdout.strip().splitlines() if l.startswith("{")][-1])
            row.update(runtime_s=m["runtime_s"], peak_mem_mb=m["peak_mem_mb"],
                       taper_size=m["taper_size"], coh_num=m["coh_num"])
            d = np.load(out_npz)
            curves[(nw_label, technique)] = d["coh"]
            freq_holder.setdefault("freq", d["freq"])
            print(f"     {m['runtime_s']}s  {m['peak_mem_mb']} MB  K={m['taper_size']}  coh_num={m['coh_num']}", flush=True)
            os.remove(out_npz)
    except Exception:
        row["error"] = traceback.format_exc()[-600:]
        print(f"     EXC: {row['error'][-200:]}", flush=True)
    metrics.append(row)
    return row


# ---- Section 1: four techniques, production configs ----
run("single-taper", NW_PROD, "s1_st", "sec1", "NW_prod(~10.8)", "plain-fft")
run("MspecBestK", NW_PROD, "s1_bk", "sec1", "NW_prod(~10.8)", "Wband=0.001")
run("FastMspec", NW_PROD, "s1_fm", "sec1", "NW_prod(~10.8)", "Wband=0.001")
run("Mspec", 100.0, "s1_ms", "sec1", "NW=100 (batch baseline)", "NW=100/K=80")

freq = freq_holder.get("freq")
if freq is not None:
    sec1_npz = {"freq": freq, "dist_km": pmd.dist_km}
    cn = next((r["coh_num"] for r in metrics if r["section"] == "sec1" and r["coh_num"]), 0)
    sec1_npz["coh_num"] = cn
    for tech in ["single-taper", "MspecBestK", "FastMspec", "Mspec"]:
        key = ("NW_prod(~10.8)", tech) if tech != "Mspec" else ("NW=100 (batch baseline)", tech)
        if key in curves:
            sec1_npz[f"coh_{SAFE[tech]}"] = curves[key]
    np.savez(OUT + "/nb3_sec1_skrh.npz", **sec1_npz)
    print("wrote nb3_sec1_skrh.npz", flush=True)

# ---- Section 1c: FastMspec vs MspecBestK cost scaling with NW ----
for nw in NW_SWEEP:
    label = f"NW={nw:g}" + (" (past ceiling)" if nw > 30 else "")
    run("FastMspec", nw, f"sw_fm_{nw:g}", "sec1c", label, f"NW={nw:g}")
    run("MspecBestK", nw, f"sw_bk_{nw:g}", "sec1c", label, f"NW={nw:g}")

sweep_npz = {"freq": freq} if freq is not None else {}
for (lbl, tech), arr in curves.items():
    if lbl.startswith("NW=") and "100" not in lbl:
        sweep_npz[f"coh_{SAFE[tech]}_{lbl.split()[0].replace('NW=', 'NW')}"] = arr
if sweep_npz:
    np.savez(OUT + "/nb3_sec1_nwsweep.npz", **sweep_npz)
    print("wrote nb3_sec1_nwsweep.npz  keys:", list(sweep_npz), flush=True)

with open(OUT + "/nb3_sec1_metrics.json", "w") as f:
    json.dump(metrics, f, indent=2)
for npy in (s1_npy, s2_npy):
    os.path.exists(npy) and os.remove(npy)
print("DONE.", flush=True)
print(json.dumps(metrics, indent=2))
