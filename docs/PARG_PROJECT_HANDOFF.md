# PARG inhibitor discovery — complete project handoff

**Investigator:** Muhammad Ali Khan
**Target:** Human poly(ADP-ribose) glycohydrolase (PARG) · UniProt Q86W56 · EC 3.2.1.143 · ChEMBL CHEMBL1795143
**Status:** analysis complete, manuscript in preparation, not submitted
**This document compiled:** 16 October 2026

---

## 0. How to use this document

This is a handoff brief for an AI assistant continuing analysis on this project. It is
self-contained: you do not have access to the underlying files, so every number you need
is reproduced here.

**Three rules for whoever picks this up:**

1. **Do not invent numbers.** Every figure in this document came from a verified run. If
   you need a quantity that is not here, say it is not here and ask for it. Several
   numbers in the project's own lab record were found to be wrong; guessing adds more.
2. **Respect the metric distinction in §5.** Some numbers are on a superseded RMSD metric
   and are explicitly flagged. Mixing them silently produces conclusions that look
   consistent and are not.
3. **Check §9 before concluding anything.** It lists what is genuinely uncertain. The
   project's main weakness is small sample size, not any single calculation.

---

## 1. The project in one paragraph

Two parallel studies on one target. **Study A** benchmarks machine-learning models for
PARG inhibitor potency prediction, screens a 3.2-million-compound commercial library with
the best model, and then runs a random control that shows the model's apparent advantage
does not survive scrutiny. **Study B** is a structural docking study showing that five
crystallographic water molecules in the PARG active site determine whether docking
recovers the correct binding pose, that they act on pose *ranking* rather than pose
*generation*, and that PARG inhibitors split into two classes which respond to those
waters in opposite directions. Study B contains the novel findings; Study A is a
careful replication of known limitations of ML in drug discovery.

---

## 2. Target selection and biology

PARG is the principal catabolic enzyme for poly(ADP-ribose), mechanistically complementary
to PARP inhibition. It is synthetically lethal in DNA-damage-response-deficient
backgrounds and sensitises cells to platinum agents and PARP inhibitors. Public chemical
space is sparse relative to PARP — fewer than 1,000 public compounds.

- Tool compound: **PDD00017273**, the first-in-class chemical probe
  (James DI et al., *ACS Chem Biol* 2016;11(11):3179–3190, doi:10.1021/acschembio.6b00609)
- Clinical precedent: **IDE161** (IDEAYA Biosciences), Phase 1, **NCT05787587**

Target chosen on four sequential criteria: biological rationale → novelty relative to
existing computational screening literature → bioactivity data sufficiency for honest
scaffold-split evaluation → structural feasibility. A ten-candidate comparison table
exists in earlier project notes but **has not been recovered** and is missing from the record.

---

## 3. STUDY A — machine learning and virtual screening

### 3.1 Dataset

Source ChEMBL target CHEMBL1795143, accessed 14 October 2026. Restricted to records with a
pChEMBL value, aggregated by canonical SMILES, duplicates averaged.

| Property | Value |
|---|---|
| Unique compounds | 870 |
| Bemis–Murcko scaffolds | 323 |
| pChEMBL range | 4.06 – 8.32 |
| Mean / SD | 6.155 / 0.764 |

**Task formulation.** At a 6.0 pChEMBL threshold the set splits 510 active / 35 inactive.
Classification on that imbalance would report high accuracy for an uninformative model, so
the task was reformulated as **regression on continuous pChEMBL before any model was
trained**. This ordering matters and should be stated in any write-up.

### 3.2 Splits

| Split | Train / Test | Role |
|---|---|---|
| Scaffold (whole Murcko families held out) | 610 / 260 | primary |
| Random | 609 / 261 | contrast |
| Temporal | 862 / 8 | abandoned, test set too small |

### 3.3 Model benchmark — all on the identical 260-compound scaffold-split test set

| Model | Representation | RMSE | R² | Seeds |
|---|---|---|---|---|
| Mean predictor | none | 0.801 | 0.000 | — |
| Tanimoto 5-NN | ECFP4 (Morgan r=2, 2048 bit) | 0.509 | 0.561 | — |
| Chemprop D-MPNN | learned graph, 318K params | 0.504 ± 0.004 | 0.569 | 3 |
| **LightGBM** | ECFP4 + 217 RDKit descriptors (2,265 features) | **0.445** | **0.664** | 5 |

LightGBM: `n_estimators=500, learning_rate=0.05, num_leaves=31`.
Random split for contrast: baseline 0.411 / R² 0.689; LightGBM 0.338 / R² 0.789.

**Finding.** Gradient boosting on classical descriptors beats the graph neural network,
and the GNN is statistically indistinguishable from nearest-neighbour similarity search
despite 318,000 learned parameters.

**Important methodological correction already applied.** Chemprop initially scored RMSE
0.467 on its *own* internally generated scaffold split (n=174), appearing to beat both the
baseline and LightGBM. Forced onto the project split (n=260) via an explicit splits file,
its error rose to **0.510** — the baseline value. The apparent advantage was an easier test
set. Seeds 0/1/2 on the matched split: 0.508, 0.506, 0.499.

**Known inconsistency, unfixed:** LightGBM used 5 seeds, Chemprop 3.

### 3.4 Library screening

Final model: LightGBM retrained on all 870 compounds, hyperparameters identical to benchmark.

Library: ZINC20 2D tranches, MW bins C–J, logP bins A–E, reactivity class A, purchasability
A and B (**in stock**, not make-on-demand). 80 tranche files →
**3,204,633 unique molecules** after deduplication.

Scored 3,204,632 (one unparseable SMILES) in 200.6 min on 8 cores, ~270 molecules/s.

| | Training set | Screened library |
|---|---|---|
| Max | 8.32 | **6.73** |
| Mean | 6.155 | **4.94** |
| Min | 4.06 | 4.11 |

The model discriminates, but predicts no library compound matching known PARG inhibitors.

*Reproduction note:* `files.docking.org` rejects the default `Python-urllib` user-agent;
a browser user-agent header is required.

### 3.5 Novelty of hits

Max Tanimoto similarity (ECFP4, 2048 bit) from each of the top 5,000 hits to the nearest
of the 870 training compounds:

| Band | n | % |
|---|---|---|
| 0.0–0.3 (very novel) | 4,428 | **88.6** |
| 0.3–0.4 | 425 | 8.5 |
| 0.4–0.6 | 121 | 2.4 |
| 0.6–1.0 | 26 | 0.5 |

The top 20 hits span predicted pChEMBL 6.73 → 6.55, a range of **0.18 units — less than
half the model's own RMSE**. Ordering within the top hits is not resolvable by the model
and must not be interpreted.

### 3.6 Random control — the decisive experiment

100 top-ranked compounds vs 100 drawn at random from the same library (seed 42, top-100
excluded from the sampling pool). Docked **interleaved** so throughput drift affects both
groups equally. Identical receptor (`rec_wet.pdbqt`), 16 Å box centred (27.20, 2.84, 0.17),
exhaustiveness 16, seed 42. Selection by model score only.

| | n | mean | median | best | below −9 |
|---|---|---|---|---|---|
| top-ranked | **92** | −8.41 | −8.42 | −10.21 | 20 (21.7%) |
| random | 100 | −8.18 | −8.17 | −10.60 | 15 (15.0%) |

**The single best-docking compound of all 192 came from the RANDOM group.**

#### Statistics — re-verified from the raw 192 affinities on 16 Oct 2026

All descriptive statistics above reproduce exactly. **The p-values in the original lab
record were computed ONE-SIDED and the record does not say so.** Both are given here:

| Test | one-sided (as recorded) | two-sided |
|---|---|---|
| Mann–Whitney, top vs random | 0.0154 | **0.0308** |
| Fisher, enrichment below −9 | 0.1536 | **0.2638** |

Sidedness confirmed by exact reproduction: `alternative='less'` → 0.0154;
`alternative='greater'` → 0.1536.

#### Sensitivity to the eight failed dockings

Eight top-ranked compounds failed to produce a pose; **zero** random compounds failed.
Failed compounds are silently excluded from summary statistics, so the comparison was
repeated with them imputed at both extremes:

| Scenario | difference | p one-sided | p two-sided | verdict |
|---|---|---|---|---|
| As run (92 v 100) | +0.227 | 0.0154 | 0.031 | holds |
| Failures at −6.42 (worst observed) | +0.067 | 0.144 | **0.284** | **fails** |
| Failures at −10.21 (best observed) | +0.370 | 0.0014 | 0.003 | holds |

**The eight failures are one molecule in eight stereoisomeric forms** — MW 459.6, 33 heavy
atoms, 10 rotatable bonds, 4 rings, logP 2.00, all identical in every computed property.
The model placed all eight in its top 100, so 8% of its best picks are a single compound.
Compounds that docked successfully averaged MW 367.6, 25.6 heavy atoms, 5.0 rotatable
bonds. Ten rotatable bonds is the likely cause of failure in 3D generation or search —
which makes the conservative imputation also the more plausible one.

**Conclusion.** The model did not demonstrably outperform random selection on this library.
Without this control the study would have reported significant enrichment and a
purchasable hit list. That conclusion would have been wrong.

---

## 4. STUDY B — structure, waters, and docking

### 4.1 Receptor preparation

Primary receptor **6OAK**, 1.70 Å, chain A, 497 ordered residues (451–961), 3 chain breaks
none within 12 Å of the site. Catalytic residues Asp737, Asn740, Glu755, Glu756, Tyr795,
Phe902 all complete. Co-crystal ligand M0M, 43 atoms / 22 heavy, centroid
(27.20, 2.84, 0.17). 360 crystallographic waters.

**Water-selection rule, fixed before any docking was run:** retain waters within **3.5 Å of
both a ligand heavy atom and a protein N or O**. Five qualify: **HOH 1211, 1232, 1307,
1317, 1402**. Two hydrogen-bond to Gln754; one contact at 2.22 Å.

Two receptors built: **dry** (no waters) and **wet** (these five).

The rule was applied unmodified to a second structure (7KG1) where it gave a different
count (4 waters). No water was chosen to produce a result. This pre-registration is the
paper's defence against the circularity objection and should be stated prominently.

**Unverified:** the Gln754 / catalytic-motif assignment has not been checked against a
primary structural source.

### 4.2 Docking protocol

```
AutoDock Vina 1.2.5
exhaustiveness 16, num_modes 9, seed 42
Open Babel 3.2.1 for PDBQT conversion, polar hydrogens retained
cross-docking box: 24 x 24 x 24 A centred (27.20, 2.84, 0.17)
screening + control box: 16 x 16 x 16 A, same centre
```

The 24 Å box was found to span the full poly(ADP-ribose) substrate groove and permit escape
to a secondary site ~6 Å below the pocket at (29.44, −0.38, −5.78). **That site is not an
artefact** — it is lined by Gln754, Tyr795, Phe738, Cys872, Gly871, Gly873, Ala874, Phe875
and is part of the genuine substrate groove. All screening used the 16 Å box.

### 4.3 Pose-accuracy metric — READ THIS BEFORE USING ANY RMSD NUMBER

Open Babel reorders ligand atoms when building the PDBQT torsion tree, so naive sequential
RMSD compares non-corresponding atoms. Measured magnitude: a dry redock reported **5.73 Å
sequentially vs 3.75 Å by optimal assignment**.

**Two metrics exist in this project:**

- **Optimal assignment** (`scipy.optimize.linear_sum_assignment` over squared interatomic
  distances, cross-element pairs forbidden). Ignores connectivity. A strict **lower bound**.
  This is the OLD metric.
- **True symmetry-corrected RMSD.** PDB Chemical Component Dictionary ideal structure matched
  onto both reference and pose as a generic graph (bond orders disregarded), minimum over all
  substructure matches, **no superposition** — compared in the crystal frame. Hydrogens
  removed explicitly. This is the CURRENT metric.

**Direct comparison across all 34 cross-docking pairs:**

| Pose class | n | mean lower-bound → true shift | max |
|---|---|---|---|
| Correct (<2 Å) | 6 | +0.010 Å | 0.06 Å |
| Incorrect | 28 | +0.930 Å | 3.42 Å |

Mann–Whitney p = 0.0002. Every shift ≥ 0, confirming the bound.

**Consequence: the old metric flattered failures specifically.** Any quantity still on the
old metric is wrong in a known direction — correct poses unaffected, incorrect poses
understated.

### 4.4 Cross-docking set

All 25 human PARG structures enumerated from UniProt Q86W56 (chain A, 448–976); 24
downloaded, 6OAK is the reference. Excluded: modified residues (MSE, CSO, CSS, CME, 6WK),
buffer components; het groups retained only at 15–50 heavy atoms; phosphorus-containing
ligands classified substrate-like and excluded.

Superposition: Kabsch on matched Cα. **490–497 atoms matched, Cα RMSD 0.13–0.37 Å** — all
structures in essentially one backbone conformation.

Result: **17 drug-like inhibitors** across three independent medicinal chemistry campaigns:

| Series | Ligands |
|---|---|
| 5LHB / 6HMx | 6WJ, 7JC, 73L, 7JB, 70J |
| 6O9x / 6OAx | M0S, M0V, M1D, M0P (×2), M0M |
| 7KGx | WD7, WDA, WDD, WDG, WDJ, WDM |

Excluded: 4 substrate-like (AR6, A8P, A1R, A3R); one surface-bound AR6 in 5A7R at 13.71 Å.
All 17 converted to PDBQT with atom counts preserved; 3–7 rotatable bonds each.

### 4.5 Cross-docking result — CURRENT, true symmetry-corrected RMSD

| Receptor | Top-1 < 2 Å | Near-native pose sampled among 9 modes | Ranked first |
|---|---|---|---|
| Dry | **0 / 17** | 6 / 17 | **0 / 6** |
| Wet | **6 / 17** | 7 / 17 | **6 / 7** |

Top-1: McNemar exact **p = 0.031**; Wilcoxon signed-rank on paired RMSD **p = 0.035**.

**Per-ligand top-1 true RMSD (Å):**

| Ligand | Dry | Wet | class |
|---|---|---|---|
| 5LHB_6WJ | 5.45 | 5.98 | displacing |
| 6HMK_7JC | 4.38 | 5.32 | displacing |
| 6HML_73L | 6.54 | 6.20 | displacing |
| 6HMM_7JB | 7.67 | 7.60 | displacing |
| 6HMN_70J | 7.54 | 8.27 | displacing |
| 6O9X_M0S | 3.74 | **0.62** | bridging |
| 6O9Y_M0V | 7.55 | 6.62 | bridging |
| 6OA0_M1D | 7.01 | 7.10 | bridging |
| 6OA1_M0P | 3.68 | **0.69** | bridging |
| 6OA3_M0M | 3.68 | 2.03 | bridging |
| 6OAL_M0P | 3.66 | **0.57** | bridging |
| 7KFP_WD7 | 2.73 | **1.94** | bridging — SEE §9.2 |
| 7KG0_WDA | 9.08 | 6.43 | bridging |
| 7KG1_WDD | 6.87 | 4.67 | bridging |
| 7KG6_WDG | 3.63 | 6.42 | bridging |
| 7KG7_WDJ | 4.59 | **0.74** | bridging |
| 7KG8_WDM | 3.36 | **0.62** | bridging |

**Note 6OA3_M0M at 2.03 Å — it fails the 2 Å criterion by 0.03 Å.** Under the old
lower-bound metric it read 1.99 and counted as a success. The headline 6/17 is one ligand
and 0.03 Å away from 7/17. Describe it as borderline, not as a failure.

### 4.6 Sampling vs ranking — the central structural finding

Computing true RMSD for **all nine Vina modes**, not only the top-ranked pose, separates
whether a near-native pose was *generated* from whether it was *selected*.

**Sampling is statistically indistinguishable.** 6/17 dry vs 7/17 wet, and the dry set is a
strict **subset** of the wet set — only 7KG7_WDJ is sampled in wet and missed in dry.
McNemar **p = 1.0**.

**Ranking is transformed.** Restricted to the 6 ligands sampled in BOTH receptors so
denominators match:

| Ligand | dry rank | wet rank |
|---|---|---|
| 6O9X_M0S | 6 | 1 |
| 6OA1_M0P | 5 | 1 |
| 6OA3_M0M | 5 | 2 |
| 6OAL_M0P | 6 | 1 |
| 7KFP_WD7 | 4 | 1 |
| 7KG8_WDM | 3 | 1 |
| **mean** | **4.83** | **1.17** |

Better in 6/6. **Wilcoxon signed-rank p = 0.031.**

**Report the rank comparison, not the ranked-first proportion.** Binarised to ranked-first
the matched comparison is 0/6 vs 5/6, McNemar **p = 0.0625 — does not clear 0.05**. The
rank test clears it because it uses the whole ranking. Also note the Wilcoxon with 6
all-same-direction pairs returns exactly 2/2⁶ = 0.031, the **minimum attainable p at that
sample size** — the test is saturated and cannot be stronger with n=6. State this.

### 4.7 Two inhibitor classes

Distance from each superposed ligand to the nearest retained bridging water:

| Class | Distance | Ligands |
|---|---|---|
| Water-bridging | 2.59 – 2.78 Å (hydrogen bond) | M0x and WD series |
| Water-displacing | 0.82 – 1.27 Å (occupies the water position) | 6WJ, 7JC, 73L, 7JB, 70J — all 5LHB/6HMx |

**Effect of water retention by class:**

| Class | n | dry | wet | mean wet−dry |
|---|---|---|---|---|
| water-bridging | 12 | 0/12 | 6/12 | −1.76 Å |
| water-displacing | 5 | 0/5 | 0/5 | +0.36 Å |

Mann–Whitney on the water effect between classes: **p = 0.014** (two-sided).

No single fixed receptor serves both classes; water retention must be decided per chemotype.

### 4.8 Rotamer switch

χ₁ torsions across all 18 structures:

| Residue | State A | State B | Separation |
|---|---|---|---|
| Tyr795 | **+66.3 ± 3.3°** (gauche⁺) | **−168.5 ± 4.0°** (trans) | ~235° |
| Phe902 | **−155.0 ± 4.4°** | **−53.1 ± 6.7°** | ~102° |

The two residues switch **in concert in all 18 structures**; no structure shows one flipped
without the other. Nine structures occupy each state.

| State A (gauche⁺) | State B (trans) |
|---|---|
| 6OAK, 6O9X, 6OA1, 6OA3, 6OAL, 7KFP, 7KG6, 7KG7, 7KG8 | 5LHB, 6HMK, 6HML, 6HMM, 6HMN, 6O9Y, 6OA0, 7KG0, 7KG1 |

Every successfully docked ligand derives from a state A structure. **This survives the RMSD
correction** — all six current successes are state A.

**B-factor objection answered:** 6OAK's Tyr795 refines at 62.4 Å², which alone would be
dismissible as disorder. But 6HMK sits at 22.1 Å² in state B and 6OA1 at 27.5 Å² in state A
— well-ordered examples occur in both clusters. The elevated B in 6OAK is a symptom of a
residue sampling two positions.

### 4.9 Ensemble docking — attempted, unsuccessful, reported as future work

Second receptor from **7KG1** (1.65 Å, state B), 496 Cα matched, four bridging waters by the
identical 3.5 Å rule, superposed into the 6OAK frame.

| Condition | Top-1 < 2 Å (OLD METRIC) |
|---|---|
| State A (6OAK) alone | 7 / 17 |
| State B (7KG1) alone | **0 / 17** |
| Ensemble, selected by Vina score | 6 / 17 |
| Oracle (perfect selection — upper bound, not a protocol) | 7 / 17 |

Score-based selection made things *worse* — it replaced a 0.62 Å pose with a 6.55 Å one on
6O9X_M0S.

Causes systematically eliminated by test: coordinate frame mismatch (6OAK files identical to
0.000 Å; rec_B vs 6OAK 0.14 Å mean, 1.09 Å max); receptor misplacement (reference ligand
makes 141 contacts <4 Å, min 2.73 Å); box coverage (site 3.49 Å from box centre); atom loss
in conversion (7745→4811 and 7562→4698, same ratio); hydrogen asymmetry (both retain polar H
only); vacated cavity wrongly emptied of water (4 waters within 5 Å of decoy vs 1 at true
site). Conclusion: receptor B is not demonstrably broken; the box was too generous and WDD
was a poor choice of test ligand (it also fails in receptor A at 4.18 Å).

---

## 5. WHAT IS CORRECTED AND WHAT IS NOT — critical

A recomputation from the lower-bound metric to true symmetry-corrected RMSD was performed
on **16 October 2026** and applied to the **cross-docking top-1 and all-pose analysis only**.

### Corrected (use these)
- §4.5 cross-docking top-1 table
- §4.6 sampling vs ranking
- §4.7 class analysis

### NOT corrected — still on the superseded lower-bound metric

| Item | Value as recorded | Status |
|---|---|---|
| Redock validation | dry 3.75 Å / wet 0.99 Å, near-native ranked 8th vs 1st | lower bound. Direction is safe: the wet value is a correct pose so its bound is tight; the dry value can only grow. The dissociation finding strengthens under correction. **Recompute anyway.** |
| Failure analysis groups | docked n=7 / failed n=10; side-chain RMSD 0.51/2.41 vs 1.57/6.55 Å; Phe902 0.63 vs 5.02 Å; Tyr795 2.41 vs 6.35 Å; "range 3.03–6.77 Å contains no structures" | **Groups are now n=6 / n=11** because 6OA3_M0M moves across. All means and the bimodal-gap claim must be recomputed. |
| Rotamer exceptions | "the one exception: 7KG6_WDG" | **Now two** — 7KG6_WDG and 6OA3_M0M. But 6OA3 fails by 0.03 Å, so describe it as a threshold artefact, not a conformational failure. |
| Ensemble docking (§4.9) | 7/17, 0/17, 6/17, 7/17 | all on the old metric, not recomputed |

### Numbers in the original lab record that are simply WRONG and have been superseded

The lab record's §9 table reads dry 2/17, wet 7/17, sampled 7/17 both, ranked first 2/7 and
7/7. **Every cell is wrong.** It also states "the same seven ligands were sampled
successfully in both receptors" — false; it is 6 dry, a strict subset of 7 wet.

---

## 6. Independent verification performed 16 October 2026

Three checks were run against primary data rather than the lab record:

1. **Coordinate frames.** Docked poses and crystal references confirmed to share one frame.
   Wet centroid offsets 0.23–0.42 Å for the six successes, 2.94–3.83 Å for the two failures,
   tracking the RMSDs. This rules out the whole RMSD table being computed across mismatched
   frames.
2. **Water geometry.** Crystal structures re-downloaded from RCSB and re-superposed
   independently (PyMOL `align`, 0.09–0.15 Å over 428–442 Cα). The bridging range reproduced
   at 2.60–2.79 Å against the record's 2.59–2.78. **One discrepancy — see §9.2.**
3. **Control statistics.** All descriptive statistics reproduced exactly from the raw 192
   affinities; p-values were found to be one-sided (§3.6).

---

## 7. Files and scripts

### Scripts (run in numerical order from `~/parg`)
```
02_curate.py          ChEMBL curation -> 870 compounds
03_split.py           scaffold / random / temporal splits
04_baseline.py        Tanimoto 5-NN baseline
05_lightgbm.py        LightGBM benchmark
06_eval_matched.py    Chemprop matched-split evaluation
07_eval_seeds.py      Chemprop seed replicates
08_inspect_pdb.py     6OAK inspection
09_verify_site.py     catalytic residue check
10_waters.py          bridging-water selection
11_prepare_receptor.py dry/wet receptors
12_redock_rmsd.py, 13_rmsd_check.py, 14_compare.py   redocking validation
15_inventory.py, 16_superpose.py, 17_prep_ligands.py cross-docking set construction
18_crossdock.py       17 ligands x 2 receptors
19_diagnose.py, 21_which_residue.py, 22_rotamer.py   failure analysis, rotamer switch
20_get_6oak.py, 23_build_stateB.py                   structure retrieval, receptor B
24-28                 ensemble attempt and diagnosis
29_train_final.py     final model, all 870 compounds
30_get_zinc.py        ZINC in-stock library retrieval
32_screen.py          score 3.2M molecules
33_novelty.py         Tanimoto similarity of hits to training set
34_control.py         random control, 100 v 100 interleaved
35_sensitivity.py     failure analysis and imputation test
--- added 15-16 October 2026 ---
36_figures.py         figures 1, 6, 7, 8 (from transcribed numbers, NOT reproducible)
38_true_rmsd.py       true symmetry-corrected RMSD, top-1
39_fig2.py            figure 3, rebuilt from CSV
40_allposes.py        true RMSD for all 9 modes -> results/crossdock_allposes.csv
41_pose_figures.py    figure 5 pose panels (PyMOL)
42_fig6.py            figure 4, sampling vs ranking
43_fig2_control.py    figure 2, random control
```

### Key data files
```
data/parg_curated.csv              870 curated compounds
data/zinc_instock.smi              3,204,633 library molecules
results/crossdock.csv              original, OLD metric
results/crossdock_truermsd.csv     CORRECTED top-1 true RMSD, 34 rows
results/crossdock_allposes.csv     all 9 modes, true RMSD
results/control.csv                192 rows: group, zinc_id, smiles, vina_kcal
results/tmp_poses/                 306 per-mode PDBs
structures/xdock/                  superposed crystal reference ligands
models/lgbm_final.pkl              final model
```

### Figures (all eight exist as of 16 Oct 2026)
```
Fig 1  model benchmark                      from transcribed numbers, NOT reproducible
Fig 2  random control distributions         reproducible, 43_fig2_control.py
Fig 3  dry vs wet cross-docking             reproducible, 39_fig2.py
Fig 4  sampling vs ranking                  reproducible, 42_fig6.py
Fig 5  docked pose panels                   reproducible, 41_pose_figures.py
Fig 6  rotamer states                       from transcribed numbers, NOT reproducible
Fig 7  side-chain separation                NOT reproducible AND on stale n=7/n=10 grouping
Fig 8  side-chain vs RMSD                   NOT reproducible AND on stale grouping
```

### Software
```
Ubuntu (WSL2) 24.04 · Python 3.11 · AutoDock Vina 1.2.5 (static linux x86_64)
Open Babel 3.2.1 · Chemprop v2
RDKit, LightGBM, scikit-learn, SciPy, NumPy — VERSIONS NOT YET RECORDED
```

---

## 8. Verified citations (PubMed, 16 October 2026)

| For | Reference | DOI |
|---|---|---|
| Vina original | Trott O, Olson AJ. *J Comput Chem* 2010;31(2):455–461 | 10.1002/jcc.21334 |
| Vina 1.2.x | Eberhardt J et al. *J Chem Inf Model* 2021;61(8):3891–3898 | 10.1021/acs.jcim.1c00203 |
| Open Babel | O'Boyle NM et al. *J Cheminform* 2011;3:33 | 10.1186/1758-2946-3-33 |
| ChEMBL | Zdrazil B et al. *Nucleic Acids Res* 2024;52(D1):D1180–D1192 | 10.1093/nar/gkad1004 |
| ZINC20 | Irwin JJ et al. *J Chem Inf Model* 2020;60(12):6065–6073 | 10.1021/acs.jcim.0c00675 |
| Murcko scaffolds | Bemis GW, Murcko MA. *J Med Chem* 1996;39(15):2887–2893 | 10.1021/jm9602928 |
| ECFP4 | Rogers D, Hahn M. *J Chem Inf Model* 2010;50(5):742–754 | 10.1021/ci100050t |
| Chemprop D-MPNN | Yang K et al. *J Chem Inf Model* 2019;59(8):3370–3388 | 10.1021/acs.jcim.9b00237 |
| PARG chemical probe | James DI et al. *ACS Chem Biol* 2016;11(11):3179–3190 | 10.1021/acschembio.6b00609 |

**Still needed, NOT verified — do not cite without checking:** LightGBM (Ke et al., NeurIPS
2017), Chemprop v2 software paper (Heid et al.), RDKit (software, cite version), Kabsch
(*Acta Cryst* 1976 + 1978 correction), MoleculeNet (Wu et al., *Chem Sci* 2018), the 2 Å
RMSD convention, NCT05787587.

---

## 9. OPEN QUESTIONS AND KNOWN WEAKNESSES

### 9.1 Sample size is the binding constraint
17 ligands. The mechanism claim rests on **6 matched pairs** with a test at its significance
floor. The water-displacing class is **n=5** with zero successes. No amount of further
computation fixes this. Any analysis proposing to strengthen the conclusions should first
ask whether n permits it.

### 9.2 UNRESOLVED: 7KFP_WD7 class assignment
Independent recomputation on 16 Oct 2026 puts 7KFP_WD7 at **0.96 Å from HOH 1402** — a
clash, i.e. water-*displacing*. The lab record assigns the entire WD series to the
water-*bridging* class. Per-water distances for 7KFP: 1211:5.40, 1232:3.25, 1307:3.44,
1317:2.11, **1402:0.96**. So WD7 appears to *displace* one water while *bridging* another,
which complicates the claim that no single receptor serves both classes.

**Impact if reclassified:**

| | bridging | displacing | Mann–Whitney p |
|---|---|---|---|
| as recorded | 0/12 → 6/12, −1.76 Å | 0/5 → 0/5, +0.36 Å | 0.0136 |
| 7KFP moved | 0/11 → 5/11, −1.85 Å | 0/6 → 1/6, +0.17 Å | 0.0145 |

The class conclusion holds either way. What is lost is the clean "0/5 displacing" line.
**This must be checked against the project's own `10_waters.py` output** — the independent
check used PyMOL `align`, the project used Kabsch on matched Cα, and that is the one place
the two could legitimately differ.

### 9.3 UNRESOLVED: test sidedness
The random control is one-sided, the cross-docking tests are two-sided, and neither is
declared. A reviewer recomputing the control will get 0.031, not 0.015. Decide and state it.

### 9.4 Multiple comparisons
At least five p-values fall between 0.01 and 0.05 across the project with no correction
applied. Nominal values are reported. A reviewer will add them up.

### 9.5 No experimental validation
All predictions are hypotheses. This caps the publication venue regardless of how clean the
analysis is.

### 9.6 Other open items
- Ten-candidate target selection table not recovered
- Gln754 catalytic-motif assignment unverified against a primary source
- LightGBM 5 seeds vs Chemprop 3 — mismatched
- 7KG6_WDG: a state A structure with no water conflict (2.59 Å) that nonetheless fails,
  carrying 6 rotatable bonds — the most of any ligand in a rigid pocket. Attributed to
  search budget. **Untested — raise exhaustiveness to confirm.**
- Figures 1, 6, 7, 8 not reproducible from CSVs; 7 and 8 also on stale grouping
- Zenodo deposition not done

---

## 10. Novelty and publication assessment

**The ML findings are confirmatory, not novel.** That deep models fail to beat similarity
baselines on molecular property prediction, that scaffold splitting degrades performance,
and that virtual screening hits fall outside the training distribution are established
results in the literature. This work replicates them cleanly on a new target.

**The novel content is PARG-specific and structural:** the five-water selection rule and the
recognition-versus-sampling dissociation; two inhibitor classes responding in opposite
directions to water retention; the concerted Tyr795/Phe902 rotamer switch.

**Realistic venues:** *J Mol Graph Model*, *Molecular Diversity*, *Comput Biol Chem*,
*J Biomol Struct Dyn*, *Scientific Reports*, MDPI journals. Out of reach: *J Med Chem*,
*JCIM*, anything selective — capped by n=17 and the absence of experimental validation.

**Framing recommendation:** the lab record's §17 currently leads with the ML negative result
as "the central finding." That is the weakest claim because it is the replication. Leading
with the PARG structural work and demoting the ML result to a supporting methodological
caution puts the novel content in front.

**The project's real strength** is methodological honesty: a random control that overturned
its own headline result and was kept in the record, plus eight documented failure modes that
each produce "wrong but publishable-looking" results. This is well above the norm for this
tier and should be foregrounded, not buried.

---

## 11. Suggested next analyses

Ordered by value per unit effort:

1. **Resolve §9.2** — recompute 7KFP_WD7's water distance with the project's own Kabsch
   superposition. Cheap, and it decides a table in the paper.
2. **Recompute the failure analysis** on the corrected n=6 / n=11 grouping, including whether
   the bimodal gap survives.
3. **Recompute the redock** (§4.9 and the 3.75/0.99 Å pair) with true RMSD.
4. **Decide and declare sidedness**; re-report the control two-sided.
5. **Test the 7KG6_WDG hypothesis** — raise exhaustiveness and see whether it is search
   budget as claimed.
6. **Regenerate figures 1, 6, 7, 8** from CSVs so the whole figure set is reproducible.
7. Record software versions; recover the target-selection table; verify Gln754; match seed
   counts; deposit on Zenodo.

Do not attempt to strengthen the statistical claims by adding tests. With n=17 and n=6 the
informative tests have already been run, and adding more worsens the multiple-comparison
problem described in §9.4.
