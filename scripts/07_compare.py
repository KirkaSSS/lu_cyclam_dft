"""Step 7: compare a new dft.jsonl with the deposited results/dft/dft.jsonl.
python 07_compare.py ../work/dft.jsonl [--ref ../results/dft/dft.jsonl] [--tol 0.05]
Matches structures by file name (a '_restrained' suffix is ignored) and reports the difference
in total energy (PBE0 + D3) in kcal/mol. Exit code 1 if any difference exceeds --tol."""
import argparse, json, os, sys
p = argparse.ArgumentParser(); p.add_argument("new"); p.add_argument("--ref", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "dft", "dft.jsonl"))
p.add_argument("--tol", type=float, default=0.05, help="kcal/mol"); a = p.parse_args()
key = lambda f: os.path.basename(f).replace("_restrained", "")
load = lambda f: {key(r["f"]): r for r in (json.loads(l) for l in open(f) if l.strip())}
ref, new = load(a.ref), load(a.new); bad = 0
print(f"{'structure':22s} {'deposited (Eh)':>18s} {'new (Eh)':>18s} {'diff (kcal/mol)':>16s}")
for k in sorted(new):
    if k not in ref: print(f"{k:22s} {'-':>18s} {new[k]['E']+new[k]['Ed3']:18.8f}   (no deposited value)"); continue
    e0, e1 = ref[k]["E"] + ref[k]["Ed3"], new[k]["E"] + new[k]["Ed3"]; d = (e1 - e0) * 627.509
    flag = "" if abs(d) <= a.tol else "  <-- exceeds tolerance"; bad += bool(flag)
    print(f"{k:22s} {e0:18.8f} {e1:18.8f} {d:16.4f}{flag}")
print("OK: all within %.2f kcal/mol" % a.tol if not bad else f"{bad} structure(s) outside tolerance")
sys.exit(1 if bad else 0)
