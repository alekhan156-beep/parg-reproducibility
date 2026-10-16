# PARG virtual screen — run package

Two-stage virtual screen against human poly(ADP-ribose) glycohydrolase (PARG).
Stage 1 scores a 3.2-million-compound library with a trained model. Stage 2 docks
the survivors into a validated receptor.

**This is CPU work, not GPU work.** Neither RDKit nor AutoDock Vina uses a graphics
card. What matters is core count. A GPU becomes useful later, for molecular
dynamics — not for anything in this package.

Stage 1 takes under an hour on 8 cores and stage 2 a few hours, so an ordinary
laptop can run the whole thing overnight. More cores simply make it faster.

---

## Files you need

From this package:
- `screen_parg.py` — stage 1, ligand-based scoring
- `dock_hits.py` — stage 2, docking

Sent separately by Ali:
- `lgbm_final.pkl` — the trained model (1.5 MB)
- `rec_wet.pdbqt` — the validated PARG receptor
- `vina_1.2.5_linux_x86_64` — the docking binary (or download it yourself)

Downloaded by you (see below):
- `zinc_instock.smi` — the compound library, ~3.2 M molecules

---

## Setup

```bash
conda create -n parg python=3.11 -y
conda activate parg
pip install rdkit lightgbm numpy
conda install -c conda-forge openbabel -y        # stage 2 only
```

Vina, if you don't have it:

```bash
wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64
chmod +x vina_1.2.5_linux_x86_64
```

---

## Getting the library

Either take `zinc_instock.smi` from Ali, or rebuild it — ZINC20 2D tranches,
molecular-weight bins C–J, logP bins A–E, reactivity class A, purchasability A and B
(**in stock**, not make-on-demand):

```python
import os, urllib.request, urllib.error
op = urllib.request.build_opener()
op.addheaders = [('User-Agent', 'Mozilla/5.0 (X11; Linux x86_64)')]
urllib.request.install_opener(op)                 # docking.org blocks Python-urllib

os.makedirs('zinc', exist_ok=True)
for m in 'CDEFGHIJ':
    for l in 'ABCDE':
        for p in 'AB':
            n = f"{m}{l}A{p}"
            if os.path.exists(f"zinc/{n}.smi"):
                continue
            try:
                urllib.request.urlretrieve(
                    f"https://files.docking.org/2D/{m}{l}/{n}.smi", f"zinc/{n}.smi")
                print(n)
            except urllib.error.HTTPError:
                pass

seen = set()
with open('zinc_instock.smi', 'w') as out:
    for f in sorted(os.listdir('zinc')):
        for line in open(f'zinc/{f}'):
            q = line.split()
            if len(q) >= 2 and q[0].lower() != 'smiles' and q[0] not in seen:
                seen.add(q[0]); out.write(f"{q[0]} {q[1]}\n")
print(len(seen), "unique molecules")
```

Expected: **3,204,633 molecules**, ~250 MB of tranche files.

---

## Stage 1 — score the library

```bash
python screen_parg.py \
  --library zinc_instock.smi \
  --model lgbm_final.pkl \
  --top 5000 \
  --jobs 16
```

Set `--jobs` to your core count. Progress prints a live rate and ETA.

**Runtime.** The model uses 2,048 fingerprint bits plus 217 RDKit descriptors.
Measured at **7.4 ms per molecule** on one core for realistic drug-like input
(the descriptors are ~99 % of that; fingerprints alone are 0.08 ms):

| Cores | Approximate wall time |
|---|---|
| 1 | ~6.6 hours |
| 8 | ~50 minutes |
| 16 | ~25 minutes |
| 32 | ~13 minutes |

Note: RDKit prints a deprecation warning per molecule unless logging is disabled.
This script disables it. If you re-enable logging, terminal I/O becomes the
bottleneck and the run appears ~25x slower than it is.

Memory stays near 1 GB whatever the library size — molecules are featurised in
chunks of 50,000 and only the scores are kept.

Outputs `screen_ranked.csv` (rank, ZINC ID, SMILES, predicted pChEMBL) and
`screen_ranked.smi` (input for stage 2).

Safe to interrupt and restart — but note stage 1 does not resume, it restarts.
Stage 2 does resume.

---

## Stage 2 — dock the survivors

```bash
python dock_hits.py \
  --hits screen_ranked.smi \
  --receptor rec_wet.pdbqt \
  --vina ./vina_1.2.5_linux_x86_64 \
  --n 500
```

Roughly 10–60 s per compound depending on cores, so 500 compounds is 2–8 hours.
Results are written after every compound, so an interrupted run loses nothing and
re-running skips what's already done.

Outputs `docking_results.csv`, sorted best-first, plus one pose file per compound in
`docked/`.

---

## Two things that will bite you if ignored

**Do not change the docking box.** It is 16 Å centred on (27.20, 2.84, 0.17). A
larger box spans the whole poly(ADP-ribose) substrate groove, and Vina will park
small ligands in a subpocket ~6 Å away that is not the inhibitor site. This was
measured, not guessed — a 24 Å box moved one test ligand 9.53 Å off; 16 Å brought it
back to 3.64 Å.

**Do not rank final hits on Vina score alone.** In validation on this exact target,
the receptor giving the *better* binding score gave the *worse* pose (−7.79 kcal/mol
at 3.75 Å versus −7.15 kcal/mol at 0.99 Å). Score and geometric correctness are
dissociated here. The ranking in `docking_results.csv` is a starting point for
inspection, not a verdict.

---

## Why the receptor has waters in it

`rec_wet.pdbqt` retains five crystallographic waters that bridge the ligand to the
protein. This is deliberate and was validated across 17 co-crystallised inhibitors:
with the waters, 7/17 docked correctly; stripped of them, 2/17. The waters do not
help Vina *find* the right pose — the same 7 were found either way — they let it
*recognise* the right pose once found. Please don't "clean up" the receptor.

---

## What to send back

- `screen_ranked.csv`
- `docking_results.csv`
- the `docked/` folder

Thanks for the compute.
