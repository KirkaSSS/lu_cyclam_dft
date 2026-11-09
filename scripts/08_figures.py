"""Step 8: Python versions of Fig. S1 (structures) and Fig. S2 (relative energies).
python 08_figures.py [--results ../results] [--out ../figures]
The figures in the paper were drawn with matlab/lu_dft_figures.m from the same files;
this script gives equivalent figures for readers without MATLAB."""
import argparse, json, os
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ase.io import read
col={71:'#2a9d8f',7:'#3b5bdb',8:'#e03131',6:'#868e96',1:'#f1f3f5'}
rad={71:0.62,7:0.30,8:0.30,6:0.27,1:0.15}
def orient(a):
    p=a.positions.copy(); Z=a.numbers; lu=list(Z).index(71); p-=p[lu]
    N=p[Z==7]; c=N.mean(0); u,s,vt=np.linalg.svd(N-c); n=vt[2]
    if np.dot(-c,n)<0: n=-n          # normal points from ring to Lu
    e1=vt[0]; e2=np.cross(n,e1)
    R=np.array([e1,n,e2])            # x=e1, y=up, z=depth
    q=p@R.T
    # small tilt to see ring
    t=np.radians(8); Rx=np.array([[1,0,0],[0,np.cos(t),-np.sin(t)],[0,np.sin(t),np.cos(t)]])
    return q@Rx.T, Z, lu
def draw(ax,a,title):
    q,Z,lu=orient(a)
    D0=np.linalg.norm(q[:,None]-q[None],axis=2)
    keep=[i for i in range(len(Z)) if not (Z[i]==1 and min(D0[i,j] for j in range(len(Z)) if Z[j]==6)<1.2)]
    q=q[keep]; Z=Z[keep]; lu=list(Z).index(71); d=np.linalg.norm(q[:,None]-q[None],axis=2); n=len(Z)
    seg=[]
    for i in range(n):
        for j in range(i+1,n):
            zi,zj=Z[i],Z[j]
            if 71 in (zi,zj):
                o=j if zi==71 else i
                if Z[o] in (7,8) and d[i,j]<3.0: seg.append((i,j,'coord'))
            elif 1 in (zi,zj):
                if d[i,j]<1.2: seg.append((i,j,'cov'))
            elif d[i,j]<1.75: seg.append((i,j,'cov'))
    # H-bond of outer-sphere water O to coordinated water H
    dl=np.linalg.norm(q-q[lu],axis=1)
    outer=[i for i in range(n) if Z[i]==8 and dl[i]>3.0]
    for o in outer:
        for h in range(n):
            if Z[h]==1 and 1.5<d[o,h]<2.3: seg.append((o,h,'hb'))
    items=[]
    for i,j,k in seg: items.append(((q[i,2]+q[j,2])/2,'b',i,j,k))
    for i in range(n): items.append((q[i,2],'a',i,None,None))
    items.sort(key=lambda x:x[0])
    for z,t,i,j,k in items:
        if t=='b':
            st=dict(cov=('-','#495057',2.2),coord=('-','#2a9d8f',2.2),hb=(':','#495057',1.6))[k]
            ax.plot(q[[i,j],0],q[[i,j],1],linestyle=st[0],color=st[1],lw=st[2],solid_capstyle='round',zorder=1+z*0.01+5)
        else:
            ax.add_patch(plt.Circle(q[i,:2],rad[Z[i]],fc=col[Z[i]],ec='#212529',lw=0.7,zorder=1+z*0.01+5.005))
    for o in outer: ax.text(q[o,0]+0.75,q[o,1]+0.55,'H$_2$O (2nd sphere)',fontsize=8.5,color='#495057',zorder=50)
    ax.set_aspect('equal'); ax.axis('off')
    ax.set_xlim(q[:,0].min()-0.8,q[:,0].max()+2.2); ax.set_ylim(-4.6,4.6)

def main():
    p = argparse.ArgumentParser(); p.add_argument("--results", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results"))
    p.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures")); a = p.parse_args()
    os.makedirs(a.out, exist_ok=True); X = os.path.join(a.results, "xtb_structures")
    R = [json.loads(l) for l in open(os.path.join(a.results, "dft", "dft.jsonl")) if l.strip()]
    C = {os.path.basename(r["f"])[:-4]: r["E"] + r["Ed3"] for r in R if r["nao"] > 400}
    e0 = min(C.values()); dE = {k: (v - e0) * 627.509 for k, v in C.items()}
    fig, axs = plt.subplots(1, 2, figsize=(10, 5.2))
    draw(axs[0], read(os.path.join(X, "s3_2_10.xyz")), "")
    draw(axs[1], read(os.path.join(X, "cn9_s4_1_6.xyz")), "")
    from matplotlib.lines import Line2D
    h = [Line2D([], [], marker="o", ls="", mfc=col[z], mec="#212529", ms=9, label=l) for z, l in [(71, "Lu"), (7, "N"), (8, "O"), (6, "C"), (1, "H (C–H omitted)")]]
    h += [Line2D([], [], color="#2a9d8f", lw=2.2, label="Lu–N/O bond"), Line2D([], [], color="#495057", ls=":", lw=1.6, label="H-bond")]
    fig.legend(handles=h, loc="lower center", ncol=7, frameon=False, fontsize=9)
    plt.tight_layout(rect=(0, 0.07, 1, 0.93))
    fig.text(0.03, 0.95, f"(a) CN 8, {dE['s3_2_10']:.1f} kcal/mol", fontsize=12)
    fig.text(0.53, 0.95, f"(b) CN 9, +{dE['cn9_s4_1_6']:.1f} kcal/mol", fontsize=12)
    plt.savefig(os.path.join(a.out, "fig1_structures_py.png"), dpi=300); plt.close(fig)
    rows = [("CN 8 (4N + 4O + H$_2$O), A", "s3_2_10"), ("CN 8 (4N + 4O + H$_2$O), B", "s4_1_6"), ("CN 9 (4N + 5O), B", "cn9_s4_1_6"),
            ("CN 9 (4N + 5O), A", "cn9_s3_2_10"), ("CN 8 (3N + 5O)", "s5_0_6")]
    fig, ax = plt.subplots(figsize=(7.2, 3.2)); y = list(range(len(rows)))[::-1]
    for (lab, k), yy in zip(rows, y):
        v = dE[k]; ax.barh(yy, max(v, 0.15), height=0.58, color="#2a6fdb" if v == 0 else "#9aa1a9")
        ax.text(v + 0.35, yy, f"{v:.1f}" if v else "0.0 (reference)", va="center", fontsize=10, color="#212529")
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=10)
    ax.set_xlabel("ΔE, PBE0-D3(BJ)/def2-SVP, ddCOSMO water (kcal/mol)", fontsize=10, color="#495057"); ax.set_xlim(0, 22)
    for s in ["top", "right"]: ax.spines[s].set_visible(False)
    for s in ["left", "bottom"]: ax.spines[s].set_color("#adb5bd")
    ax.tick_params(colors="#495057", length=0); plt.tight_layout()
    plt.savefig(os.path.join(a.out, "fig2_energies_py.png"), dpi=300); plt.close(fig)
    print("Figures written to", os.path.abspath(a.out))
if __name__ == "__main__":
    main()
