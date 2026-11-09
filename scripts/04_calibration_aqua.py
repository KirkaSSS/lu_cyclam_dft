"""Korak 4: kalibracija na [Lu(H2O)9]3+ (CN 9) naspram [Lu(H2O)8]3+ + H2O (CN 8+1).  python 04_calibration_aqua.py [--out ../work]"""
from lib import *

def aq(nin,nout,seed):
    rng=np.random.default_rng(seed); pos=[np.zeros(3)];Z=[71]
    d=rng.normal(size=(nin,3)); d/=np.linalg.norm(d,axis=1)[:,None]
    for _ in range(3000):
        g=-(d[:,None]-d[None])/(np.linalg.norm(d[:,None]-d[None],axis=2)[...,None]+np.eye(nin)[...,None])**3
        d=d-0.01*(g*(1-np.eye(nin))[...,None]).sum(1); d/=np.linalg.norm(d,axis=1)[:,None]
    Os=[]
    def w(o,u):
        u=u/np.linalg.norm(u);a=np.cross(u,rng.normal(size=3));a/=np.linalg.norm(a)
        return [o,o+.97*(.62*u+.79*a),o+.97*(.62*u-.79*a)]
    for u in d:
        Os.append(2.35*u); 
        for q,z in zip(w(2.35*u,u),[8,1,1]):pos.append(q);Z.append(z)
    for j in range(nout):
        for _ in range(999):
            o=Os[rng.integers(nin)]; u=rng.normal(size=3);u/=np.linalg.norm(u);o=o+2.7*u
            if np.linalg.norm(o)>4 and min(np.linalg.norm(o-p) for p in pos)>2.2:break
        for q,z in zip(w(o,o),[8,1,1]):pos.append(q);Z.append(z)
    return Atoms(numbers=Z,positions=np.array(pos))
import argparse, os
from ase.io import write
p = argparse.ArgumentParser(); p.add_argument("--out", default="../work"); a = p.parse_args(); os.makedirs(a.out, exist_ok=True)
for nin, nout, s in [(9, 0, 1), (8, 1, 0)]:
    m = aq(nin, nout, s); m.calc = X(); LBFGS(m, logfile=None, maxstep=.08).run(fmax=.03, steps=500)
    d = m.get_all_distances()[0]; O = sorted(round(float(d[i]), 2) for i, z in enumerate(m.numbers) if z == 8)
    write(os.path.join(a.out, f"aq{nin}_{nout}.xyz"), m); print(nin, nout, round(m.get_potential_energy() / Hartree, 5), O, flush=True)
