"""Korak 1: GFN2-xTB (ALPB voda) pretraga konformera [Lu(cyclam)(H2O)n]3+ + (5-n) H2O.
Upotreba:  python 01_conformers.py [--seeds 16] [--jobs 2] [--out ../work]"""
import argparse, json, os, time
from concurrent.futures import ProcessPoolExecutor
def one(args):
    nin, nout, seed, out = args
    from lib import build, X, LBFGS, Hartree, analyse
    a = build(nin, nout, seed); a.calc = X()
    t = time.time(); opt = LBFGS(a, logfile=None, maxstep=0.1); opt.run(fmax=0.03, steps=600)
    E = a.get_potential_energy() / Hartree; N, O, co, cn = analyse(a)
    name = f"s{nin}_{nout}_{seed}"; a.write(os.path.join(out, name + ".xyz"))
    return dict(name=name, start=f"{nin}+{nout}", E_h=float(E), steps=int(opt.get_number_of_steps()),
                CN_N=int(cn), CN_O=int(co), LuN=[float(x) for x in N], LuO=[float(x) for x in O], t=round(time.time() - t))
if __name__ == "__main__":
    p = argparse.ArgumentParser(); p.add_argument("--seeds", type=int, default=16)
    p.add_argument("--jobs", type=int, default=2); p.add_argument("--out", default="../work")
    a = p.parse_args(); os.makedirs(a.out, exist_ok=True)
    tasks = [(i, o, s, a.out) for s in range(1, a.seeds + 1) for i, o in [(5, 0), (4, 1), (3, 2)]]
    with ProcessPoolExecutor(a.jobs) as ex, open(os.path.join(a.out, "xtb_conformers.jsonl"), "w") as f:
        for r in ex.map(one, tasks):
            f.write(json.dumps(r) + "\n"); f.flush(); print(r["name"], r["CN_N"], r["CN_O"], round(r["E_h"], 5), flush=True)
