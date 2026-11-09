"""Ceo tok od pocetka: python run_all.py [--seeds 16] [--jobs 2] [--threads 2] [--skip-xtb]
--skip-xtb preskace korake 1-4 i koristi geometrije iz ../results/xtb_structures (iste kao u analizi)."""
import argparse, os, subprocess, sys
p = argparse.ArgumentParser(); p.add_argument("--seeds", type=int, default=16); p.add_argument("--jobs", type=int, default=2)
p.add_argument("--threads", type=int, default=2); p.add_argument("--memory", type=int, default=5500)
p.add_argument("--skip-xtb", action="store_true"); a = p.parse_args()
here = os.path.dirname(os.path.abspath(__file__)); os.chdir(here)
work = os.path.abspath("../work"); os.makedirs(work, exist_ok=True)
env = dict(os.environ, OMP_NUM_THREADS="1")
run = lambda *c: subprocess.run([sys.executable, *c], check=True, env=env)
if a.skip_xtb:
    geo = os.path.abspath("../results/xtb_structures"); suffix = ""
else:
    run("01_conformers.py", "--seeds", str(a.seeds), "--jobs", str(a.jobs), "--out", work)
    run("02_rank_xtb.py", "--dir", work)
    run("03_cn9_cages.py", "--dir", work)
    run("04_calibration_aqua.py", "--out", work)
    geo = work; suffix = "_restrained"
# napomena: s novim seed-ovima najbolje strukture mogu imati druga imena -> proverite xtb_ranking.json
targets = ["aq9_0", "aq8_1", f"cn9_s4_1_6{suffix}", "s4_1_6", "s3_2_10", "s5_0_6", f"cn9_s3_2_10{suffix}"]
out = os.path.join(work, "dft.jsonl"); env["OMP_NUM_THREADS"] = str(a.threads)
for t in targets:
    run("05_dft_single_point.py", os.path.join(geo, t + ".xyz"), "--threads", str(a.threads), "--memory", str(a.memory), "--out", out)
run("06_summary.py", out)
