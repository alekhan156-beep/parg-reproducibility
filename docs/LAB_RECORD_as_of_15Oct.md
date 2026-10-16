# PARG inhibitor discovery — project record

**Investigator:** Muhammad Ali Khan
**Target:** Human poly(ADP-ribose) glycohydrolase (PARG)
**UniProt:** Q86W56 · **EC:** 3.2.1.143 · **ChEMBL target:** CHEMBL1795143
**Working directory:** `~/parg/`
**Record compiled:** 15 October 2026

This file is the authoritative record of what was done, in what order, with what
settings, and what came out. Every number here was produced on the investigator's
machine and appears in terminal output. Manuscript text should be written from
this file, not from memory.

---

## 0. Division of labour

Docking execution and molecular dynamics are performed by the investigator.
Analysis scripting, record-keeping and manuscript drafting are collaborative.

---

## 1. Software environment

| Component | Version | Role |
|---|---|---|
| Ubuntu (WSL2) | 24.04 | operating environment |
| Python | 3.11 | all analysis |
| AutoDock Vina | 1.2.5 (static linux x86_64) | docking |
| Open Babel | 3.2.1 | PDBQT conversion |
| RDKit | *record before submission* | descriptors, fingerprints |
| LightGBM | *record before submission* | gradient boosting |
| Chemprop | v2 | D-MPNN |
| scikit-learn / SciPy / NumPy | *record before submission* | baselines, assignment, arrays |

**Action required:** run the version dump and fill the blanks above.

---

## 2. Target selection

Criteria applied in sequence: biological rationale → novelty (not already saturated
by computational screening literature) → bioactivity data sufficiency (enough
compounds across enough scaffolds for honest scaffold-split evaluation) →
structural feasibility (high-resolution co-crystals, complete catalytic machinery,
multiple independent complexes for protocol validation).

PARG satisfied all four. **The ten-candidate comparison table is not yet in this
record — it exists in earlier project notes and must be recovered and inserted.**

Supporting context for the introduction:
- PARG is the principal catabolic enzyme for poly(ADP-ribose); mechanistically
  complementary to PARP inhibition
- Synthetically lethal in DNA-damage-response-deficient backgrounds
- Sensitises to platinum agents and PARP inhibitors
- Tool compound: PDD00017273 (first-in-class chemical probe)
- Clinical precedent: **IDE161 (IDEAYA Biosciences), Phase 1, NCT05787587**
- Chemical space sparse relative to PARP — fewer than 1,000 public compounds

---

## 3. Dataset

Source: ChEMBL target CHEMBL1795143, accessed 14 October 2026.
Records restricted to those with a pChEMBL value, aggregated by canonical SMILES,
duplicates averaged.

| Property | Value |
|---|---|
| Unique compounds | 870 |
| Bemis–Murcko scaffolds | 323 |
| pChEMBL range | 4.06 – 8.32 |
| Mean / SD | 6.155 / 0.764 |
| File | `data/parg_curated.csv` |

### Decision: regression, not classification
At a 6.0 pChEMBL threshold the set partitions 510 active / 35 inactive. Classification
metrics on that imbalance would report high accuracy for an uninformative model.
Task reformulated as regression on continuous pChEMBL **before any model was trained**.
Script `02_curate.py` was edited accordingly (`agg['label'] = agg['pchembl']`).

### Splits (`03_split.py`)
| Split | Train / Test | Status |
|---|---|---|
| Scaffold (whole Murcko families held out) | 610 / 260 | **primary** |
| Random | 609 / 261 | contrast |
| Temporal | 862 / 8 | abandoned — test set too small |

---

## 4. Model benchmark — all on the identical 260-compound scaffold-split test set

| Model | Representation | RMSE | R² | Seeds |
|---|---|---|---|---|
| Mean predictor | none | 0.801 | 0.000 | — |
| Tanimoto 5-NN | ECFP4 (Morgan r=2, 2048 bit) | 0.509 | 0.561 | — |
| Chemprop D-MPNN | learned molecular graph, 318K params | 0.504 ± 0.004 | 0.569 | 3 |
| **LightGBM** | ECFP4 + 217 RDKit descriptors (2265 features) | **0.445** | **0.664** | 5 |

Random split for contrast: baseline 0.411 / R² 0.689; LightGBM 0.338 / R² 0.789.
The scaffold-minus-random gap is the cost of extrapolating to unseen chemotypes.

LightGBM settings: `n_estimators=500, learning_rate=0.05, num_leaves=31`.

### Headline finding
Gradient boosting on classical descriptors beats the graph neural network, and the
graph neural network is statistically indistinguishable from nearest-neighbour
similarity search despite 318K learned parameters.

### Correction made — matched test sets
Chemprop initially scored **0.467** on its own internal scaffold split (n=174),
appearing to beat both the baseline and LightGBM. Forced onto the project split
(n=260) via an explicit splits file its error rose to **0.510** — the baseline value.
The apparent advantage was an easier test set. Seeds 0/1/2 on the matched split:
0.508, 0.506, 0.499.

Chemprop splits-file format (v2) is a list of dicts:
`[{'train': [...], 'val': [...], 'test': [...]}]`

**Known inconsistency to fix:** LightGBM used 5 seeds, Chemprop 3. Make them match.

---

## 5. Structure preparation

Primary receptor: **6OAK**, 1.70 Å.

| Property | Value |
|---|---|
| Chain | A |
| Ordered residues | 497 (451–961) |
| Chain breaks | 3, none within 12 Å of the site |
| Catalytic residues present | Asp737, Asn740, Glu755, Glu756, Tyr795, Phe902 — all complete |
| Co-crystal ligand | M0M, 43 atoms / 22 heavy |
| Ligand centroid | (27.20, 2.84, 0.17) |
| Crystallographic waters | 360 |

### Water selection rule — applied before any docking
Retain waters within **3.5 Å of both a ligand heavy atom and a protein N or O**.
Five qualify: **HOH 1211, 1232, 1307, 1317, 1402**. Two hydrogen-bond to Gln754,
the conserved glutamine preceding the catalytic motif; one contact at 2.22 Å.

Two receptors built: dry (no waters) and wet (these five).

**Action required:** verify the Gln754 / catalytic-motif assignment against a primary
structural paper before it appears in the manuscript.

---

## 6. Docking protocol

```
AutoDock Vina 1.2.5
box 24 × 24 × 24 Å centred on (27.20, 2.84, 0.17)
exhaustiveness 16, num_modes 9, seed 42
```

### Pose accuracy metric — important methodological note
Open Babel reorders ligand atoms when building the PDBQT torsion tree, so naïve
sequential RMSD compares non-corresponding atoms. All RMSD values in this record
are computed by **optimal assignment under an element-identity constraint**
(`scipy.optimize.linear_sum_assignment` over squared interatomic distances,
cross-element pairs forbidden).

This metric is a **lower bound** on true symmetry-corrected RMSD — it ignores
connectivity and can only underestimate. Thresholds should be read conservatively.

**Action required:** recompute with RDKit `GetBestRMS` before submission.

Magnitude of the artefact: dry-receptor redock reported **5.73 Å** sequentially
versus **3.75 Å** by assignment.

---

## 7. Redocking validation (native M0M into 6OAK)

| Receptor | Best affinity | Top-1 RMSD | Rank of near-native pose |
|---|---|---|---|
| Dry | −7.79 kcal/mol | 3.75 Å | **8 of 9** (1.15 Å) |
| Wet (5 waters) | −7.15 kcal/mol | **0.99 Å** (centroid offset 0.48 Å) | **1** |

### Finding: score and pose accuracy are dissociated
The receptor giving the *better* score gave the *worse* pose. Binding score is not a
proxy for geometric correctness. Virtual screening output cannot be ranked on raw
score without independent pose validation.

---

## 8. Cross-docking set construction

All 25 human PARG structures enumerated from UniProt Q86W56 (all chain A, residues
448–976). 24 downloaded (6OAK is the reference).

Filtering: modified residues (MSE, CSO, CSS, CME, 6WK) and buffer components
excluded; het groups retained only at 15–50 heavy atoms; phosphorus-containing
ligands classified substrate-like (ADP-ribose analogues) and excluded from the
inhibitor analysis.

Superposition: Kabsch on matched Cα atoms. **490–497 atoms matched, Cα RMSD
0.13–0.37 Å** — all structures in essentially the same backbone conformation.

Result: **17 drug-like inhibitors in the active site**, spanning three independent
medicinal chemistry campaigns:

| Series | Ligands |
|---|---|
| 5LHB / 6HMx | 6WJ, 7JC, 73L, 7JB, 70J |
| 6O9x / 6OAx | M0S, M0V, M1D, M0P (×2), M0M |
| 7KGx | WD7, WDA, WDD, WDG, WDJ, WDM |

Excluded: 4 substrate-like (AR6, A8P, A1R, A3R); one surface-bound AR6 in 5A7R at
13.71 Å from the site.

All 17 converted to PDBQT with atom counts preserved; 3–7 rotatable bonds each.

---

## 9. Cross-docking result (24 Å box)

| Receptor | Top-1 < 2 Å | Near-native pose sampled | Sampled pose ranked first |
|---|---|---|---|
| Dry | **2 / 17** | 7 / 17 | 2 / 7 |
| Wet | **7 / 17** | 7 / 17 | **7 / 7** |

The **same seven ligands** were sampled successfully in both receptors.

### Finding: bridging waters govern recognition, not sampling
Water retention did not change which poses Vina could generate — only whether it
could identify them. In the wet receptor every near-native pose generated was ranked
first; in the dry receptor only two of seven were.

Median top-1 is **not** a valid summary here — the distribution is bimodal.

### Per-ligand top-1 RMSD (Å)

| Ligand | Dry | Wet |
|---|---|---|
| 5LHB_6WJ | 4.98 | 5.25 |
| 6HMK_7JC | 3.88 | 4.81 |
| 6HML_73L | 5.50 | 5.67 |
| 6HMM_7JB | 6.10 | 5.35 |
| 6HMN_70J | 5.08 | 7.82 |
| 6O9X_M0S | 2.56 | **0.62** |
| 6O9Y_M0V | 7.04 | 6.47 |
| 6OA0_M1D | 5.62 | 4.71 |
| 6OA1_M0P | 2.58 | **0.69** |
| 6OA3_M0M | **1.91** | **1.99** |
| 6OAL_M0P | 3.46 | **0.57** |
| 7KFP_WD7 | **1.96** | **1.88** |
| 7KG0_WDA | 5.66 | 6.35 |
| 7KG1_WDD | 6.26 | 4.18 |
| 7KG6_WDG | 3.55 | 5.98 |
| 7KG7_WDJ | 3.93 | **0.74** |
| 7KG8_WDM | 3.11 | **0.62** |

---

## 10. Failure analysis

Binding-site residues defined as those within 6 Å of M0M (23 residues).

| Group | Mean side-chain RMSD | Mean max displacement |
|---|---|---|
| Docked < 2 Å (n=7) | 0.51 Å | 2.41 Å |
| Failed (n=10) | 1.57 Å | 6.55 Å |

**The range 3.03 – 6.77 Å contains no structures.** Complete bimodal separation.

### Residues responsible

| Residue | Movement, successes | Movement, failures | Gap |
|---|---|---|---|
| **Phe902** (catalytic) | 0.63 Å | 5.02 Å | 4.39 Å |
| **Tyr795** (catalytic) | 2.41 Å | 6.35 Å | 3.94 Å |
| Cys872 | 0.25 | 1.16 | 0.91 |
| all others | — | — | ≤ 0.91 Å |

---

## 11. Rotamer switch — the principal structural finding

χ₁ torsions measured across all 18 structures:

| Residue | State A | State B | Separation |
|---|---|---|---|
| Tyr795 | **+66.3 ± 3.3°** (gauche⁺) | **−168.5 ± 4.0°** (trans) | ~235° |
| Phe902 | **−155.0 ± 4.4°** | **−53.1 ± 6.7°** | ~102° |

- The two residues switch **in concert in all 18 structures**. No structure shows one
  flipped without the other.
- Nine structures occupy each state.
- Every ligand successfully docked derives from a state A structure. Every
  conformationally-driven failure derives from state B.

### B-factor objection answered
6OAK's Tyr795 refines at 62.4 Å², which alone would be dismissible as disorder. But
**6HMK sits at 22.1 Å² in state B and 6OA1 at 27.5 Å² in state A** — well-ordered
examples occur in both clusters. The elevated B in 6OAK is a symptom of a residue
sampling two positions, not an argument against the finding.

### Structure assignment

| State A (gauche⁺) | State B (trans) |
|---|---|
| 6OAK, 6O9X, 6OA1, 6OA3, 6OAL, 7KFP, 7KG6, 7KG7, 7KG8 | 5LHB, 6HMK, 6HML, 6HMM, 6HMN, 6O9Y, 6OA0, 7KG0, 7KG1 |

### The one exception
**7KG6_WDG** is a state A structure (Tyr795 +66.7°, Phe902 −146.9°) with no water
conflict (2.59 Å) that nonetheless failed to dock. It carries six rotatable bonds, the
most of any ligand in a rigid pocket. Attributed to search budget, not conformation.
**Untested — raise exhaustiveness to confirm.**

---

## 12. Second structural finding — two inhibitor classes

Distance from each superposed ligand to the nearest retained bridging water:

| Class | Distance | Ligands |
|---|---|---|
| Water-bridging | 2.59 – 2.78 Å (hydrogen-bond distance) | M0x and WD series |
| **Water-displacing** | **0.82 – 1.27 Å** (occupying the water position) | 6WJ, 7JC, 73L, 7JB, 70J — all 5LHB/6HMx |

Two distinct binding strategies at one site. **No single fixed receptor can serve
both classes**; water retention must be decided per chemotype, not globally.

---

## 13. Ensemble docking attempt — unsuccessful, reported as future work

Second receptor built from **7KG1** (1.65 Å, state B), 496 Cα matched, four bridging
waters by the identical 3.5 Å rule, superposed into the 6OAK frame.

| Condition | Top-1 < 2 Å |
|---|---|
| State A (6OAK) alone | 7 / 17 |
| State B (7KG1) alone | **0 / 17** |
| Ensemble, selected by Vina score | 6 / 17 |
| Oracle (perfect selection — upper bound, not a protocol) | 7 / 17 |

Score-based selection made things *worse* — it replaced a 0.62 Å pose with a 6.55 Å
one on 6O9X_M0S.

### Causes eliminated by test
| Hypothesis | Test result |
|---|---|
| Coordinate frames mismatched | 6OAK files identical to 0.000 Å; rec_B vs 6OAK 0.14 Å mean, 1.09 Å max |
| Receptor misplaced | reference ligand makes 141 contacts < 4 Å, minimum 2.73 Å |
| Box doesn't cover the site | site 3.49 Å from box centre |
| Atoms lost in conversion | 7745→4811 and 7562→4698, same ratio |
| Hydrogen asymmetry between receptors | both deposited with H, both retain polar H only |
| Vacated cavity wrongly emptied of water | 4 waters within 5 Å of decoy vs 1 at true site — not a water-filled void |

### What was established
- The 24 Å box lets Vina escape to a site ~6 Å below the pocket, at
  (29.44, −0.38, −5.78). Tightening to 16 Å moved WDD from 9.53 Å to 3.64 Å offset.
- **That site is not an artefact.** It is lined by Gln754, Tyr795, Phe738, Cys872,
  Gly871, Gly873, Ala874, Phe875 — i.e. part of the genuine poly(ADP-ribose)
  substrate groove, which extends well beyond where a 22-atom inhibitor sits.
- **WDD was a poor choice of self-redock test ligand** — it also fails in receptor A
  (4.18 Å, 3.49 Å centroid offset), so it could not isolate a receptor-B problem.

### Conclusion
Receptor B is not demonstrably broken; the box was too generous and the test ligand
poorly chosen. A 16 Å re-run was started and abandoned when throughput collapsed
(1 of 34 jobs in 20 minutes, likely thermal throttling). **Ensemble treatment is
reported as indicated future work.**

**Carried forward as mandatory:** reduce the box to 16 Å before library screening.
A 24 Å box spans the whole substrate groove and will contaminate screening output at
scale with compounds that score well in the wrong subpocket.

---

## 14. Library screening

### 14.1 Final model
Retrained on all 870 compounds (no held-out set; benchmarking complete), identical
hyperparameters to the benchmark. 2,265 features. Saved as `models/lgbm_final.pkl`.

### 14.2 Library
ZINC20 2D tranches, `https://files.docking.org/2D/{MW}{logP}/`. Molecular-weight
bins C–J, logP bins A–E, reactivity class A (anodyne), purchasability A and B
(**in stock** — orderable, not make-on-demand). 80 tranche files, all retrieved.

**3,204,633 unique molecules** after deduplication → `data/zinc_instock.smi`.

*Technical note:* docking.org rejects the default `Python-urllib` user-agent. A
browser user-agent header is required or every request returns an HTTP error.

### 14.3 Scoring
3,204,632 molecules scored (one unparseable SMILES) in 200.6 min on 8 cores,
~270 molecules/s.

| | Training set | Screened library |
|---|---|---|
| Max | 8.32 | **6.73** |
| Mean | 6.155 | **4.94** |
| Min | 4.06 | 4.11 |

The model discriminates — a 1.8-unit spread is over two training standard
deviations — but predicts no library compound to match known PARG inhibitors, and
the whole distribution sits below the training mean. Given that the training set is
the product of dedicated medicinal chemistry against this target, a general catalogue
containing better compounds would have been the surprising result.

---

## 15. Novelty of screening hits

Maximum Tanimoto similarity (ECFP4, 2048 bit) from each of the top 5,000 hits to the
nearest of the 870 training compounds:

| Band | n | % |
|---|---|---|
| Tanimoto 0.0–0.3 (very novel) | 4,428 | **88.6** |
| 0.3–0.4 (novel) | 425 | 8.5 |
| 0.4–0.6 (related) | 121 | 2.4 |
| 0.6–1.0 (close analogue) | 26 | 0.5 |

The top 20 hits all fall between 0.16 and 0.30 — chemically unrelated to anything the
model trained on.

**Consequence.** The model's measured accuracy (RMSE 0.445) was established on
PARG-like chemistry under scaffold splitting. These predictions are far outside that
regime; the true error is larger and unquantified.

**A direct indication of this:** the top 20 hits span 6.73 to 6.55, a range of
**0.18 units — less than half the model's own RMSE**. The ordering within the top
hits is not resolvable by the model and should not be interpreted.

---

## 16. Random control — the decisive experiment

### 16.1 Design
100 top-ranked compounds against 100 drawn at random from the same library (seed 42,
top-100 excluded from the sampling pool). Docked **interleaved**, so that machine
throttling over the run affects both groups equally. Identical receptor
(`rec_wet.pdbqt`), box (16 Å, centred 27.20 / 2.84 / 0.17), exhaustiveness 16, seed 42.

Selection by model score only; no RMSD or outcome information entered the choice.

### 16.2 Result as run

| | n | mean | median | best | below −9 |
|---|---|---|---|---|---|
| top-ranked | **92** | −8.41 | −8.42 | −10.21 | 20 (21.7%) |
| random | 100 | −8.18 | −8.17 | −10.60 | 15 (15.0%) |

Difference +0.227 kcal/mol, Mann–Whitney U **p = 0.015**.

Note that the single best-docking compound of all 192 came from the **random** group.

### 16.3 The eight failures
Eight top-ranked compounds failed to produce a docked pose; **zero** random compounds
failed. Their properties:

| ZINC ID | MW | heavy | rotB | rings | logP |
|---|---|---|---|---|---|
| 107487267, 107487263, 107487261, 107487269 | 459.6 | 33 | 10 | 4 | 2.00 |
| 102402661, 102402652, 102402654, 102402657 | 459.6 | 33 | 10 | 4 | 2.00 |

**All eight are identical in every computed property.** They are one molecule in eight
stereoisomeric forms. The model placed all eight in its top 100 — 8% of its best picks
are a single compound.

Against the compounds that docked successfully (mean MW 367.6, 25.6 heavy atoms,
5.0 rotatable bonds), the failures are substantially larger and twice as flexible.
Ten rotatable bonds is the likely cause of failure in 3D generation or search.

### 16.4 Sensitivity analysis
Because failed compounds are silently excluded, the comparison was repeated with the
eight imputed at the extremes of the observed range:

| Scenario | difference | p | |
|---|---|---|---|
| As run (92 v 100) | +0.227 | 0.0154 | holds |
| Failures imputed at −6.42 (worst observed) | +0.067 | **0.144** | **fails** |
| Failures imputed at −10.21 (best observed) | +0.370 | 0.0014 | holds |

**The result does not survive the conservative assumption.** The conservative
assumption is also the more plausible one: those compounds failed *because* they are
large and highly flexible, which is precisely the profile that docks poorly.

### 16.5 Enrichment
20/92 (21.7%) of top-ranked versus 15/100 (15.0%) of random compounds scored below
−9 kcal/mol: **1.45-fold enrichment, Fisher exact p = 0.154 — not significant.**

### 16.6 Conclusion
**The model did not demonstrably outperform random selection on this library.** The
nominally significant Mann–Whitney result rests on the exclusion of eight compounds
and does not survive conservative imputation; the enrichment metric was never
significant. This is a null result, not a demonstration of failure — the evidence is
insufficient to claim an advantage.

Without this control the study would have reported significant enrichment
(p = 0.015) and a purchasable hit list. That conclusion would have been wrong.

---

## 17. Synthesis — the central finding

Four independent measurements in this study point the same direction:

1. **Chemprop D-MPNN (0.504 ± 0.004) is statistically indistinguishable from
   Tanimoto nearest-neighbour search (0.509)** despite 318K learned parameters.
2. **Scaffold-split error substantially exceeds random-split error** for every model
   (LightGBM 0.445 vs 0.338; baseline 0.509 vs 0.411), quantifying the cost of
   unfamiliar chemistry.
3. **88.6% of screening hits fall below Tanimoto 0.3** to any training compound, and
   the model cannot resolve differences among its own top hits.
4. **No demonstrable advantage over random selection** when ranking a library outside
   the training distribution.

Taken together: these models operate by recognising chemistry resembling their
training data, and their advantage degrades toward zero as that resemblance falls
away. This is a measurable, bounded claim supported by four separate lines of
evidence within one study.

---

## 18. Methodological controls — these belong in the manuscript

Each is a failure mode that produces a wrong but publishable-looking result.

1. **Test-set matching.** Models compared on independently generated splits are not
   comparable. Demonstrated: 0.467 → 0.510 for the same model.
2. **Atom ordering in pose comparison.** Sequential RMSD is invalid after PDBQT
   conversion. Demonstrated: 5.73 → 3.75 Å.
3. **Task formulation.** 510/35 class imbalance makes classification uninformative.
4. **Water selection precedes docking.** The geometric rule was fixed before any
   docking ran and applied unmodified to a second structure, where it gave a different
   count (4 vs 5). No water was chosen to produce a result.
5. **Ensemble selection by score, not RMSD.** Selection by RMSD would use the answer
   being tested. Oracle figures reported separately and labelled as upper bounds.
6. **Score ≠ accuracy.** Demonstrated twice: better affinity, worse pose.
7. **Random control on the virtual screen.** Model-ranked compounds compared against
   randomly selected ones under identical conditions, interleaved to control for
   machine drift. **This control overturned the headline screening result.**
8. **Sensitivity analysis for dropped compounds.** Failed dockings are silently
   excluded from summary statistics; imputing them at the extremes tests whether a
   conclusion depends on which compounds happened to fail. Here it did.

---

## 19. Limitations

- RMSD reported as a lower bound (see §6)
- Protonation by Open Babel; no pKa prediction, no tautomer enumeration
- 17 inhibitors from 3 campaigns — narrow chemical sampling, several close analogues
- Receptor flexibility treated as two discrete rigid states; no continuous relaxation,
  no backbone motion, no induced fit on binding
- Random control used 100 per group; a larger control would narrow the confidence
  interval on the null result
- Screening used a single model (LightGBM); Random Forest and XGBoost were available
  in the environment but not benchmarked, and neither was an uncertainty estimate
  derived for individual predictions
- No experimental validation — all predictions are hypotheses

---

## 20. Outstanding actions before submission

- [ ] Recover and insert the ten-candidate target selection table
- [ ] Record remaining software versions (environment.yml and requirements.txt
      captured from the `parg` environment — verify RDKit, LightGBM, Chemprop entries)
- [ ] Recompute RMSD with RDKit `GetBestRMS`
- [ ] Verify the Gln754 / catalytic motif claim against a primary source
- [ ] Match seed counts between LightGBM (5) and Chemprop (3)
- [ ] Regenerate all figures from CSVs by a script in `scripts/` (currently built from
      transcribed numbers — not reproducible)
- [ ] Add a figure for the random control: score distributions, both groups
- [ ] Consider benchmarking Random Forest for per-prediction uncertainty
- [ ] Deposit scripts and data on Zenodo; cite the DOI in Methods

---

## 21. Script index

| Script | Purpose |
|---|---|
| 02_curate.py | ChEMBL curation → 870 compounds |
| 03_split.py | scaffold / random / temporal splits |
| 04_baseline.py | Tanimoto 5-NN baseline |
| 05_lightgbm.py | LightGBM benchmark |
| 06_eval_matched.py, 07_eval_seeds.py | Chemprop matched-split evaluation |
| 08_inspect_pdb.py, 09_verify_site.py | 6OAK inspection, catalytic residue check |
| 10_waters.py, 11_prepare_receptor.py | bridging-water selection, dry/wet receptors |
| 12_redock_rmsd.py, 13_rmsd_check.py, 14_compare.py | redocking validation |
| 15_inventory.py, 16_superpose.py, 17_prep_ligands.py | cross-docking set construction |
| 18_crossdock.py | 17 ligands × 2 receptors |
| 19_diagnose.py, 21_which_residue.py, 22_rotamer.py | failure analysis, rotamer switch |
| 20_get_6oak.py, 23_build_stateB.py | structure retrieval, receptor B |
| 24–28 | ensemble attempt and diagnosis |
| 29_train_final.py | final model, all 870 compounds |
| 30_get_zinc.py | ZINC in-stock library retrieval |
| 32_screen.py | score 3.2M molecules |
| 33_novelty.py | Tanimoto similarity of hits to training set |
| 34_control.py | random control, 100 v 100 interleaved |
| 35_sensitivity.py | failure analysis and imputation test |

Scripts run in numerical order reproduce every number in this record.

---

*Record updated 15 October 2026, 23:00 GST, after completion of the random control.*
