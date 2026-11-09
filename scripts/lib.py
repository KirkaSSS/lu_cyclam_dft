import numpy as np, sys
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.optimize import LBFGS
from ase.units import Bohr, Hartree
from tblite.interface import Calculator as TB
from rdkit import Chem
from rdkit.Chem import AllChem

class X(Calculator):
    implemented_properties=['energy','forces']
    def __init__(s,charge=3,solv='water',**k):
        super().__init__(**k); s.ch=charge; s.sv=solv
    def calculate(s,atoms=None,properties=['energy'],system_changes=all_changes):
        super().calculate(atoms,properties,system_changes)
        c=TB("GFN2-xTB",atoms.numbers,atoms.positions/Bohr,charge=s.ch,uhf=0)
        c.set("verbosity",0)
        if s.sv: c.add("alpb-solvation",s.sv)
        r=c.singlepoint()
        s.results['energy']=r.get('energy')*Hartree
        s.results['forces']=-r.get('gradient')*Hartree/Bohr

def cyclam(seed):
    m=Chem.AddHs(Chem.MolFromSmiles('C1CNCCNCCCNCCNC1'))
    AllChem.EmbedMolecule(m,randomSeed=seed)
    AllChem.MMFFOptimizeMolecule(m)
    return m

def build(nin,nout,seed):
    """Lu on the N4 axis; nin waters capping, nout H-bonded second-sphere. Ring relaxed from a metal-bound guess."""
    rng=np.random.default_rng(seed)
    m=cyclam(seed); p=m.GetConformer().GetPositions(); Z=[a.GetAtomicNum() for a in m.GetAtoms()]
    N=[i for i,z in enumerate(Z) if z==7]; cen=p[N].mean(0)
    u,sv,vt=np.linalg.svd(p[N]-cen); nrm=vt[2]
    p=p-cen; lu=nrm*1.5
    pos=list(p); sym=list(Z)+[71]; pos.append(lu)
    # capping waters: Fibonacci-like points on the cap opposite N4
    def water(o,d):
        d=d/np.linalg.norm(d)
        a=np.cross(d,rng.normal(size=3)); a/=np.linalg.norm(a)
        c,s_=np.cos(np.radians(52)),np.sin(np.radians(52))
        return [o,o+0.97*(c*d+s_*a),o+0.97*(c*d-s_*a)]
    dirs=[]
    if nin:
        ang=rng.uniform(0,6.28)
        if nin>=5: dirs=[nrm]+[np.array(0) for _ in range(0)]
        t=np.radians(65)
        ring=nin-1 if nin>=4 else nin
        e1=np.cross(nrm,[1,0.2,0.3]); e1/=np.linalg.norm(e1); e2=np.cross(nrm,e1)
        dirs=[nrm*np.cos(t)+np.sin(t)*(np.cos(ang+2*np.pi*i/ring)*e1+np.sin(ang+2*np.pi*i/ring)*e2) for i in range(ring)]
        if nin>=4: dirs.append(nrm)
    Os=[]
    for d in dirs:
        o=lu+2.35*d; Os.append(o)
        for q,z in zip(water(o,d),[8,1,1]): pos.append(q); sym.append(z)
    for j in range(nout):
        # H-bond acceptor at ~2.7 A from a random inner O, outward, rejecting clashes
        for _ in range(500):
            b=Os[rng.integers(len(Os))]; d=rng.normal(size=3); d/=np.linalg.norm(d)
            o=b+2.7*d
            if np.dot(o-lu,nrm)>0 and np.linalg.norm(o-lu)>4.0 and min(np.linalg.norm(o-np.array(x)) for x in pos)>2.2:
                break
        for q,z in zip(water(o,o-lu),[8,1,1]): pos.append(q); sym.append(z)
    return Atoms(numbers=sym,positions=np.array(pos))

def analyse(a,f=sys.stdout):
    lu=[i for i,z in enumerate(a.numbers) if z==71][0]
    d=a.get_all_distances()[lu]
    N=sorted(round(d[i],2) for i,z in enumerate(a.numbers) if z==7)
    O=sorted(round(d[i],2) for i,z in enumerate(a.numbers) if z==8)
    return N,O,sum(x<3.0 for x in O),sum(x<3.0 for x in N)
