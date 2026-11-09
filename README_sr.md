# Koordinacija Lu(III) u Lu–ciklam kompleksu: GFN2-xTB / DFT podaci i kod

**Autor:** Snežana Đurković · ORCID [0009-0007-6638-0682](https://orcid.org/0009-0007-6638-0682) · Institut za nuklearne nauke „Vinča“, Univerzitet u Beogradu

Paket sadrži kod, strukture i rezultate nezavisne DFT provere iz Supporting Information
(Section S1) rada o Lu(III) kompleksu sa ciklamom (1,4,8,11-tetraazaciklotetradekan).
Proverava se da li je katjon [Lu(cyclam)(H2O)5]3+ nonakoordinovan (4N + 5O) ili
oktakoordinovan (4N + 4O), sa petim molekulom vode u drugoj koordinacionoj sferi.
Provera dopunjuje ωB97X-D4/def2-TZVP optimizacije iz glavnog teksta.

Glavni rezultat: najniža je oktakoordinovana struktura [Lu(cyclam)(H2O)4]3+·H2O; najbolji
nonakoordinovan izomer je 6,9 kcal/mol iznad nje (PBE0-D3(BJ)/def2-SVP, ddCOSMO voda,
na GFN2-xTB geometrijama).

**Nonakoordinovan izomer:** struktura istog sastava, [Lu(cyclam)(H2O)5]3+, u kojoj je Lu(III) vezan za
devet donorskih atoma: četiri atoma azota iz ciklama i atome kiseonika svih pet molekula vode
(koordinacioni broj 9; 4N + 5O). Izomer je oktakoordinovanog oblika jer oba sadrže iste atome;
razlikuju se samo po položaju jednog molekula vode: vezan je za Lu ili je vodoničnom vezom vezan u
drugoj koordinacionoj sferi (4N + 4O + H2O). U kodu se donorski atom smatra vezanim ako mu je
rastojanje od Lu manje od 3,0 Å.

## Sadržaj

| Putanja | Sadržaj |
| --- | --- |
| `results/xtb_structures/` | Sve GFN2-xTB strukture (.xyz, Å; u drugom redu je GFN2-xTB energija u eV) |
| `results/dft/dft.jsonl` | DFT rezultati, jedan red po strukturi (energije u hartree) |
| `results/verification/` | Izlazi provere reproduktivnosti (vidi dole) |
| `figures/` | Objavljene Fig. S1 i S2 (iz `matlab/lu_dft_figures.m`), PDF i PNG |
| `scripts/` | Python koraci 01–08 (vidi dole) |
| `matlab/` | MATLAB funkcije: tabela rezultata i objavljene slike |
| `requirements.txt` | Tačne verzije korišćenih paketa |
| `README.md` | Uputstvo na engleskom |
| `SHA256SUMS.txt` | Kontrolne sume svih fajlova |

## Model i ulazni podaci

Račun ne koristi nijedan eksperimentalni podatak. Ulaz je hemijski model katjona: jedan Lu3+,
jedan ciklam (SMILES u `scripts/lib.py`) i pet molekula vode, ukupno naelektrisanje +3.
Kontrajoni (BF4−, Na+) nisu uključeni; voda je modelovana kao kontinuum (ε = 78,36).
Sve početne geometrije pravi kod; nasumični brojevi koriste fiksne početne vrednosti (seed).

## Postupak

| Korak | Skripta | Šta radi |
| --- | --- | --- |
| 1 | `01_conformers.py` | 16 konformera ciklama (RDKit, MMFF94) × 3 raspodele vode (5+0, 4+1, 3+2 vezane + spolja) = 48 struktura; GFN2-xTB optimizacija u ALPB vodi (f_max < 0,03 eV/Å) |
| 2 | `02_rank_xtb.py` | Rangira sve strukture po GFN2-xTB energiji (`xtb_ranking.json`) |
| 3 | `03_cn9_cages.py` | Prisilni 4N + 5O kavezi od najboljih struktura, pa relaksacija bez ograničenja |
| 4 | `04_calibration_aqua.py` | Referentni sistem: [Lu(H2O)9]3+ naspram [Lu(H2O)8]3+·H2O |
| 5 | `05_dft_single_point.py` | PBE0/def2-SVP, ECP28MWB na Lu, density fitting, ddCOSMO voda, + D3(BJ); Mulliken naboj na Lu i dipolni momenti |
| 6 | `06_summary.py` | Tabela relativnih energija |
| 7 | `07_compare.py` | Poredi novi `dft.jsonl` sa deponovanim |
| 8 | `08_figures.py` | Python verzija Fig. S1–S2 (objavljene slike su iz `matlab/lu_dft_figures.m`) |

`run_all.py` pokreće korake 1–6 redom.

## Zahtevi

- Linux ili macOS (na Windowsu: WSL2). PySCF ne radi direktno na Windowsu.
- Python 3.10–3.13 (rezultati su dobijeni sa Python 3.13 na Linux x86_64), najmanje 8 GB RAM-a.

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

## Reprodukcija rezultata

Sve komande iz foldera `scripts/`.

**A. Pregled deponovanih rezultata (nekoliko sekundi)**

```bash
python 06_summary.py ../results/dft/dft.jsonl
python 08_figures.py
```

**B. Ponovni DFT račun na deponovanim geometrijama (~1 h 50 min na 2 jezgra)**

```bash
python run_all.py --skip-xtb --threads 4
python 07_compare.py ../work/dft.jsonl
```

**C. Ceo tok, uključujući pretragu konformera (~2 h na 2 jezgra)**

```bash
python run_all.py --threads 4 --jobs 4
python 07_compare.py ../work/dft.jsonl
```

`07_compare.py` ispisuje razliku od deponovanih energija u kcal/mol i završava se kodom 1
ako neka razlika pređe toleranciju (podrazumevano 0,05 kcal/mol).

## Očekivani rezultati

| Struktura (fajl) | Koordinacija Lu | ΔE (kcal/mol) | ΔE bez D3 (kcal/mol) |
| --- | --- | --- | --- |
| `s3_2_10` | 4N + 4O + H2O (druga sfera) | 0,0 | 0,0 |
| `s4_1_6` | 4N + 4O + H2O (druga sfera) | 4,2 | 4,2 |
| `cn9_s4_1_6` | 4N + 5O | 6,9 | 7,8 |
| `cn9_s3_2_10` | 4N + 5O | 14,9 | 15,9 |
| `s5_0_6` | 3N + 5O (jedan N odvojen) | 18,9 | 19,4 |
| `aq8_1` (referenca) | 8 H2O + 1 H2O (druga sfera) | 0,0 | 0,0 |
| `aq9_0` (referenca) | 9 H2O | 1,2 | 1,8 |

Imena fajlova: `sN_M_K` = početak sa N vezanih i M spoljašnjih molekula vode, konformer K
(konačna koordinacija može se razlikovati od početne); `cn9_*` = 4N + 5O kavez; `aq*` = Lu3+ u vodi.

## Provera reproduktivnosti

Paket je 9. 11. 2026. instaliran u novo Python 3.13 okruženje iz `requirements.txt` i ponovo
pokrenut (Linux x86_64, 2 jezgra):

- Koraci 1–4 (48 konformera, CN 9 kavezi, referentni sistem): sve GFN2-xTB energije su identične
  deponovanim strukturama; `s3_2_10` je ponovo najniža (`results/verification/xtb_ranking.json`).
- Korak 5 za `aq8_1`, `aq9_0` i ponovo dobijeni `s3_2_10`: ukupne energije su iste kao deponovane,
  do 10⁻⁸ Eh (`results/verification/compare_output.txt`).

## Napomene o reproduktivnosti

- Deponovani `cn9_*.xyz` fajlovi su geometrije sa ograničenjima, korišćene za DFT. Novi račun
  koraka 3 piše ih kao `cn9_*_restrained.xyz` (isto značenje), a uz njih i oslobođene geometrije
  `cn9_*.xyz`; oslobađanje menja GFN2-xTB energiju za manje od 0,1 kcal/mol i zadržava CN 9.
- Početne strukture koriste fiksne seed vrednosti, pa sa verzijama iz `requirements.txt` pretraga
  daje iste strukture. Druge verzije ili platforme mogu malo da promene numeriku, pa time i koji je
  konformer najniži; tada proverite `work/xtb_ranking.json` pre koraka 5. Put B od toga ne zavisi.
- Energije su elektronske, iz single-point računa na GFN2-xTB geometrijama sa double-zeta bazom,
  bez termičkih i entropijskih korekcija. GFN2-xTB služi samo za geometrije: u testu sa Lu3+ u vodi
  pogrešno favorizuje CN 9 (za 8,8 kcal/mol), dok DFT nivo daje eksperimentalni CN 8.
- Dipolni moment katjona zavisi od izbora koordinatnog početka; navedeni su oko centra naboja
  jezgara i oko Lu.

## MATLAB

Folder `matlab/` ne traži toolbox-ove (za izvoz slika: MATLAB R2020a ili noviji).

```matlab
cd lu_cyclam_dft/matlab
T = lu_dft_results;      % tabela deponovanih rezultata
lu_dft_figures           % Fig. S1 i S2 (PNG 300 dpi i PDF) u ../figures
lu_dft_setup             % jednom: Python okruženje (Windows: u WSL2)
T = lu_dft_run('SkipXtb', true, 'Threads', 4);   % put B iz MATLAB-a
```

Na Windowsu `lu_dft_setup` i `lu_dft_run` traže WSL2 sa Ubuntu-om. Instalacija, jednom:
u PowerShell-u kao administrator `wsl --install -d Ubuntu`, restart, pa u Ubuntu-u
`sudo apt update && sudo apt install -y python3-venv python3-pip`. `lu_dft_results` i
`lu_dft_figures` rade bez Pythona i WSL-a.

## Metode i programi za citiranje

GFN2-xTB: Bannwarth, Ehlert, Grimme, J. Chem. Theory Comput. 15 (2019) 1652 ·
ALPB: Ehlert et al., J. Chem. Theory Comput. 17 (2021) 4250 ·
PBE0: Adamo, Barone, J. Chem. Phys. 110 (1999) 6158 ·
D3 / D3(BJ): Grimme et al., J. Chem. Phys. 132 (2010) 154104; J. Comput. Chem. 32 (2011) 1456 ·
def2: Weigend, Ahlrichs, Phys. Chem. Chem. Phys. 7 (2005) 3297; lantanidi: Gulde, Pollak,
Weigend, J. Chem. Theory Comput. 8 (2012) 4062 · ECP28MWB: Dolg, Stoll, Preuss, J. Chem. Phys. 90
(1989) 1730 · ddCOSMO: Lipparini et al., J. Chem. Theory Comput. 9 (2013) 3637 ·
PySCF: Sun et al., J. Chem. Phys. 153 (2020) 024109 · tblite, simple-dftd3, ASE, RDKit,
Basis Set Exchange (verzije u `requirements.txt`).
