"""Korak 5: PBE0/def2-SVP (ECP28MWB na Lu), density fitting, ddCOSMO voda + D3(BJ), PySCF.
python 05_dft_single_point.py struktura.xyz [--charge 3] [--basis def2-svp] [--threads 2] [--memory 5500] [--out rezultati.jsonl]
Jedan racun za kompleks (~517 baznih funkcija) traje ~20 min na 2 jezgra i trazi ~4-5 GB RAM."""
import argparse, json, time, numpy as np
from pyscf import gto, dft, lib
from ase.io import read
import basis_set_exchange as bse
from dftd3.interface import RationalDampingParam, DispersionModel
p = argparse.ArgumentParser(); p.add_argument("xyz"); p.add_argument("--charge", type=int, default=3)
p.add_argument("--basis", default="def2-svp"); p.add_argument("--threads", type=int, default=2)
p.add_argument("--memory", type=int, default=5500, help="MB"); p.add_argument("--out", default=None)
a = p.parse_args(); lib.num_threads(a.threads)
m = read(a.xyz)
nw = bse.get_basis(a.basis, elements=[71], fmt="nwchem", header=False)   # Lu baza + ECP (def2 -> ECP28MWB)
bs = {"default": a.basis, "Lu": gto.basis.parse(nw.split("ECP")[0].replace('BASIS "ao basis" SPHERICAL PRINT', "").replace("END", ""))}
ecp = gto.basis.parse_ecp("ECP" + nw.split("ECP", 1)[1])
mol = gto.M(atom=[(s, tuple(q)) for s, q in zip(m.get_chemical_symbols(), m.positions)], basis=bs, ecp={"Lu": ecp},
            charge=a.charge, spin=0, verbose=0, max_memory=a.memory)
mf = dft.RKS(mol, xc="pbe0").density_fit(); mf.grids.level = 2; mf.conv_tol = 1e-7
ms = mf.ddCOSMO(); ms.with_solvent.eps = 78.36
t = time.time(); e = ms.kernel()
ed = DispersionModel(m.numbers, m.positions / 0.529177249).get_dispersion(RationalDampingParam(method="pbe0"), grad=False)["energy"]
dm = ms.make_rdm1(); ao = mol.aoslice_by_atom(); S = mol.intor("int1e_ovlp"); pop = (dm @ S).diagonal()
lu = m.get_chemical_symbols().index("Lu"); qM = mol.atom_charge(lu) - pop[ao[lu, 2]:ao[lu, 3]].sum()
r = mol.atom_coords(); Zq = mol.atom_charges()
def dip(o):   # dipol katjona zavisi od pocetka: racuna se oko centra naboja jezgara i oko Lu
    with mol.with_common_orig(o): i = mol.intor_symmetric("int1e_r")
    return ((Zq[:, None] * (r - o)).sum(0) - np.einsum("xij,ji", i, dm)) * 2.541746
cm = (Zq[:, None] * r).sum(0) / Zq.sum()
res = dict(f=a.xyz, E=float(e), Ed3=float(ed), E_total=float(e + ed), conv=bool(ms.converged), qLuMull=float(qM),
           dipCoC=float(np.linalg.norm(dip(cm))), dipLu=float(np.linalg.norm(dip(r[lu]))), nao=mol.nao, t=round(time.time() - t))
print(json.dumps(res), flush=True)
if a.out: open(a.out, "a").write(json.dumps(res) + "\n")
