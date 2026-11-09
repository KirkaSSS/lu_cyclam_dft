"""Korak 3: prisilni CN 9 kavezi (4 N + 5 O) sa blagim ogranicenjima, pa relaksacija bez njih.
Pise cn9_<ime>_restrained.xyz (geometrija koriscena za DFT u analizi) i cn9_<ime>.xyz (posle oslobadjanja).
python 03_cn9_cages.py [--dir ../work] [--start s5_0_3 s4_1_6 s3_2_10 s4_1_9]"""
import argparse, os
from ase.io import read, write
from ase.constraints import Hookean
from lib import X, LBFGS, Hartree, analyse
p = argparse.ArgumentParser(); p.add_argument("--dir", default="../work")
p.add_argument("--start", nargs="+", default=["s5_0_3", "s4_1_6", "s3_2_10", "s4_1_9"]); a = p.parse_args()
for f in a.start:
    m = read(os.path.join(a.dir, f + ".xyz")); m.calc = X(); lu = list(m.numbers).index(71)
    cs = [Hookean(lu, int(i), k=15., rt=2.62) for i, z in enumerate(m.numbers) if z == 8]
    cs += [Hookean(lu, int(i), k=15., rt=2.68) for i, z in enumerate(m.numbers) if z == 7]
    m.set_constraint(cs); LBFGS(m, logfile=None, maxstep=0.08).run(fmax=0.03, steps=500); m.set_constraint()
    e = m.get_potential_energy() / Hartree; N, O, co, cn = analyse(m)
    write(os.path.join(a.dir, f"cn9_{f}_restrained.xyz"), m)
    print(f, "restrained", round(e, 5), [float(x) for x in N], [float(x) for x in O], flush=True)
    b = m.copy(); b.calc = X(); LBFGS(b, logfile=None, maxstep=0.08).run(fmax=0.03, steps=400)
    N, O, co, cn = analyse(b); write(os.path.join(a.dir, f"cn9_{f}.xyz"), b)
    print(" released", round(b.get_potential_energy() / Hartree, 5), [float(x) for x in N], [float(x) for x in O], "CN", int(cn) + int(co), flush=True)
