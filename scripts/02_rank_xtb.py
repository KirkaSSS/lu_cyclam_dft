"""Korak 2: rangiranje svih GFN2-xTB struktura po energiji.  python 02_rank_xtb.py [--dir ../work]"""
import argparse, glob, json, os
from ase.io import read
from lib import X, Hartree, analyse
p = argparse.ArgumentParser(); p.add_argument("--dir", default="../work"); a = p.parse_args()
R = []
for f in glob.glob(os.path.join(a.dir, "s*_*_*.xyz")):
    m = read(f); m.calc = X(); e = m.get_potential_energy() / Hartree; N, O, co, cn = analyse(m)
    R.append(dict(file=os.path.basename(f), E_h=float(e), CN_N=int(cn), CN_O=int(co), LuN=[float(x) for x in N], LuO=[float(x) for x in O]))
R.sort(key=lambda r: r["E_h"]); e0 = R[0]["E_h"]
for r in R: print(f"{r['file']:14s} {(r['E_h']-e0)*627.509:6.1f} kcal/mol  N{r['CN_N']} O{r['CN_O']}  Lu-N {r['LuN']}  Lu-O {r['LuO']}")
json.dump(R, open(os.path.join(a.dir, "xtb_ranking.json"), "w"), indent=1)
