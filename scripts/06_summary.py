"""Korak 6: tabela relativnih DFT energija.  python 06_summary.py [rezultati.jsonl]"""
import json, sys
f = sys.argv[1] if len(sys.argv) > 1 else "../work/dft.jsonl"
R = [json.loads(l) for l in open(f) if l.strip()]
for grp, lab in [(lambda r: r["nao"] < 400, "Lu-akva kalibracija"), (lambda r: r["nao"] >= 400, "Lu-ciklam kompleks")]:
    c = [r for r in R if grp(r)]
    if not c: continue
    tot = lambda r: r["E"] + r["Ed3"]; e0 = min(map(tot, c)); e1 = min(r["E"] for r in c)
    print(f"\n{lab}\n{'struktura':26s} {'dE D3':>7s} {'dE bez D3':>9s} {'q(Lu)':>6s} {'mu CoC/Lu (D)':>14s}")
    for r in sorted(c, key=tot):
        print(f"{r['f'].split('/')[-1]:26s} {(tot(r)-e0)*627.509:7.1f} {(r['E']-e1)*627.509:9.1f} {r['qLuMull']:6.2f} {r['dipCoC']:6.1f} / {r['dipLu']:.1f}  conv={r['conv']}")
