# Computational methodology for PARG ligand modelling and docking

Muhammad Ali Khan

Methods and reproducibility record dated 16 October 2026

## 1 Study design and scope

We investigated human poly(ADP-ribose) glycohydrolase (PARG; UniProt Q86W56; ChEMBL CHEMBL1795143) using two complementary computational analyses. The first evaluated molecular representations for activity regression and used the selected model to prioritize compounds from a purchasable compound library. The second examined the effects of retained crystallographic waters and receptor conformation on ligand pose prediction. Model accuracy, pose recovery, pose ranking, and docking-score enrichment were treated as distinct outcomes. No biochemical inhibition, cellular activity, selectivity, or therapeutic efficacy was measured.

The completed historical analyses used the files in /home/aliwhoami/parg. A fresh execution was launched on 16 October 2026 in /home/aliwhoami/parg/runs/20261016_132012. This execution uses frozen data and prepared structure snapshots and produces new models, scores, docking outputs, and logs. Historical completed procedures are described below in the past tense. The separate execution appendix identifies procedures implemented for the current run and records its status at document creation; queued computations are not represented as completed experiments.

## 2 Bioactivity collection and chemical curation

### 2.1 Source records and eligibility

PARG-associated activity records were obtained from ChEMBL [1]. The laboratory record gives an access date of 14 October 2026. The preserved raw table, data/parg_raw.csv, contained 1,244 records: 666 IC50, 577 EC50, and one Kd measurement. The current execution starts from this table rather than querying the live database. An exact ChEMBL release identifier and the original retrieval command are not preserved in the inspected analysis scripts, so the frozen file and its checksum identify the reproducible input.

The curation script retained rows with a missing data_validity_comment, standard_relation equal to '=', nonmissing canonical_smiles, and nonmissing pchembl_value. The activity field was converted to a numeric value with invalid values coerced to missing, and remaining missing values were removed. These criteria retained 1,218 measurements, comprising 650 IC50, 567 EC50, and one Kd record. The code does not add an endpoint-specific, assay-format, or confidence-score restriction at this curation stage. Consequently, the regression target is pooled ChEMBL activity rather than a homogeneous biochemical pIC50 endpoint.

### 2.2 Structure processing and aggregation

SMILES strings were parsed with RDKit. The largest fragment was selected using LargestFragmentChooser, the selected fragment was processed with Uncharger, and structures containing more than 60 heavy atoms were excluded. Structures that could not be parsed or standardized were rejected. Canonical SMILES and InChIKeys were calculated for the resulting molecules. All 1,218 eligible activity records survived the recorded structure-standardization step.

Measurements were grouped by the standardized structure InChIKey. The target for each compound was the median pChEMBL value across its retained measurements. Within-compound standard deviation and record count were also retained. The first observed standardized SMILES, activity type, and assay identifier were stored, together with the earliest document year. The first assay identifier and endpoint therefore do not describe every measurement contributing to a compound median. Earlier descriptions of averaging duplicate canonical SMILES are superseded by this code-level account of InChIKey grouping and median aggregation.

The resulting dataset contained 870 compounds. The continuous median pChEMBL value was assigned to the label field for regression. Descriptive threshold counts were 510 compounds with pChEMBL at least 6, 35 with pChEMBL at most 5, and 325 between 5 and 6. These are three groups, not a two-class partition. An obsolete print statement in the original curation script compares the continuous label with 0 and 1 and prints zero counts; it does not alter the saved regression targets.

### 2.3 Scaffold construction and dataset partitions

RDKit MurckoScaffoldSmiles was used with includeChirality=False to generate scaffolds, producing 323 distinct scaffolds. For the primary partition, compounds were grouped by scaffold and groups were sorted by decreasing size. Complete groups were assigned to training while the accumulated training count remained below 70% of all compounds; subsequent groups were assigned to testing. A group was not split to meet an exact target size. This deterministic, size-ordered assignment produced 610 training and 260 test compounds. The fresh execution explicitly checked that the two sets shared no scaffold strings and reproduced the saved compound identities, targets, and assignments.

For comparison, a random partition was generated using numpy.random.default_rng(42), a permutation of all compound indices, and a cut at int(0.7 × n). This produced 609 training and 261 test compounds. A temporal partition assigned compounds with earliest document year at or before 2019 to training and all others to testing; it produced 862 and eight compounds and was not used for the principal benchmark. The implemented comparison assigns missing years to the test side. Neither the scaffold partition nor the random partition constitutes repeated cross-validation.

## 3 Molecular representations and regression models

### 3.1 Mean and similarity baselines

The mean baseline predicted the mean training-set activity for every test compound. Its R² was calculated against the variance of the test targets; it was not fixed to zero. For the similarity baseline, molecules were represented by 2,048-bit Morgan fingerprints with radius 2, corresponding to an ECFP4-style representation [2]. Chirality encoding was not explicitly enabled. The five training compounds with the largest Tanimoto similarities to each test compound were selected using NumPy sorting. Their activity values were averaged with Tanimoto similarities as weights. If all five weights summed to zero, the training mean was used.

For bit vectors A and B, Tanimoto similarity was the number of shared set bits divided by the number of bits set in either fingerprint. The five-neighbour prediction was the sum of similarity-weighted activities divided by the sum of the five similarities. The test molecules were never included in the training-neighbour pool. The historical baseline additionally printed Pearson correlation, whereas RMSE and R² were the common comparison measures.

### 3.2 LightGBM representation and fitting

For LightGBM [3], the Morgan fingerprint was concatenated with the 217 RDKit descriptors available in the installed descriptor collection, producing 2,265 input features. The original benchmark calculated descriptors with CalcMolDescriptors; the final model and library-screening scripts iterated over Descriptors._descList. NaN and positive or negative infinite values were replaced with zero. No feature scaling, feature selection, or independently documented hyperparameter search was applied.

LGBMRegressor was fitted with 500 estimators, learning rate 0.05, and 31 leaves. Five fits were performed with random_state values 0 through 4 for each partition, with other modelling parameters left at the installed package defaults. All five fits reproduced identical predictions within each partition in the fresh execution. The seed standard deviation therefore describes numerical repeatability under these settings, not uncertainty from alternative training samples or scaffold partitions. The fresh execution explicitly uses four training threads and stores the ordered descriptor names with the final model.

After benchmarking, a final LightGBM model with the same hyperparameters and random_state=0 was fitted to all 870 compounds for library scoring. This model has no held-out evaluation set of its own; the previously held-out benchmark evaluates the modelling procedure. The final model was serialized with pickle, together with feature metadata. In the new execution, the ordered descriptor list is checked before screening, preventing silent use of a changed descriptor ordering with the same feature count.

### 3.3 Directed message passing neural network

Chemprop 2.3.1 was used for regression with a directed bond-message-passing architecture [4,5]. An earlier experiment used Chemprop's own scaffold partition with 174 test compounds and is retained only as exploratory history. The comparable benchmark used data/fixed_split.json and the same 260 test compounds as the scaffold benchmark. The remaining compounds were allocated to 518 training and 92 validation examples. Thus, the test population was matched, but the number of examples contributing to gradient updates differed from the 610 examples used by LightGBM and the similarity baseline. The supplied split file was reused directly; its original train-versus-validation construction procedure is not encoded in the inspected split-generation script.

The saved configuration specifies the V2 multihot atom featurizer, message hidden dimension 300, message-passing depth 3, norm aggregation with normalization constant 100, ReLU activation, and zero dropout. The feed-forward predictor used hidden dimension 300 and one hidden layer. The model contained approximately 318,000 trainable parameters. No additional molecular descriptor file, pretrained encoder, or ensemble averaging was specified. Each run used batch size 64, zero data-loader workers, 50 epochs, and one model.

The learning-rate schedule used two warm-up epochs, initial learning rate 0.0001, maximum learning rate 0.001, and final learning rate 0.0001. Validation loss was the checkpoint-tracking metric, with min_delta=0.0. PyTorch seeds 0, 1, and 2 were evaluated independently; data_seed was 0 and num_replicates was 1. Accelerator and device selection were configured as automatic. During the new execution CUDA was unavailable, so training ran on the CPU. Training, validation, and test index sets were checked for overlap and complete coverage, and the Chemprop input SMILES and targets were checked against the curated table.

### 3.4 Evaluation and interpretation

For n test compounds, RMSE = sqrt[sum((y − prediction)²) / n] and R² = 1 − sum((y − prediction)²) / sum((y − mean(y_test))²). Models were evaluated on their held-out predictions, with means and population standard deviations across the available model seeds where applicable. Chemprop predictions were joined to the test targets by SMILES with a one-to-one merge and a complete-test-set count check in the new execution. The historical evaluation scripts merged against the full truth table without these explicit coverage checks.

The scaffold test was the principal evaluation and the random test supplied a comparison for chemical-series overlap. No confidence interval, bootstrap over scaffolds, nested validation, or formal paired model-comparison test was performed in the inspected benchmark. Numerical proximity between two RMSE values is therefore not evidence of statistical equivalence, and different seed counts do not create comparable uncertainty estimates across model families.

## 4 Library construction and virtual screening

### 4.1 ZINC library snapshot

The library was assembled from ZINC20 two-dimensional tranche files [6]. The retrieval script enumerated molecular-weight bin letters C through J, logP bin letters A through E, a fixed third character A, and final purchasability letters A and B, yielding 80 tranche names. Files were requested from https://files.docking.org/2D/{MW}{logP}/{tranche}.smi with a browser-style user agent. The laboratory record interprets the requested classes as anodyne, in-stock compounds; availability refers to the historical download and was not reconfirmed with suppliers.

Downloaded files were read in sorted filename order. Header lines and malformed rows were skipped, and exact input SMILES strings were deduplicated using a Python set. The first retained identifier for each input SMILES was preserved. This was string-level deduplication, not canonical-structure, stereoisomer-family, or InChIKey deduplication. The merged input contained 3,204,633 entries. Unlike the training curation, the screening scripts did not repeat largest-fragment selection and neutralization on every library molecule; this processing difference is part of the implemented workflow.

### 4.2 Feature calculation and ranking

Library molecules were parsed with RDKit and represented with the same 2,048 fingerprint bits and 217 descriptor functions used by the final LightGBM model. Descriptor NaN and infinite values were replaced by zero and screening features were converted to float32. The original script used eight worker processes and batches of 50,000 while retaining the full input and prediction lists in memory. It saved the 5,000 highest predictions to a CSV and SMILES file. The historical record reports 3,204,632 scored entries and one unsuccessful entry; the original function returns the same missing result for a parsing failure and a descriptor exception, so the exact failure cause is not independently distinguished by that output.

The fresh run uses the same fingerprint and descriptor calculation with eight workers, batches of 4,096, and chunksize 64. It streams all successful scores to screen_all.csv, writes rejected entries to screen_failures.csv, and retains a heap of the top 5,000. Exact score ties are resolved by earlier input order. This tie convention is explicit and may differ from the original NumPy argsort ordering for equal scores. Progress is written after each batch. Every saved prediction is a model estimate of the pooled activity target rather than a measured affinity.

### 4.3 Similarity to known compounds

For each of the 5,000 highest-ranked entries, the maximum Morgan-fingerprint Tanimoto similarity to the 870 curated compounds was calculated. Descriptive intervals were [0,0.3), [0.3,0.4), [0.4,0.6), and [0.6,1.01). These intervals summarize similarity and do not constitute a calibrated applicability-domain model. Low similarity was not taken as proof of patent novelty, new mechanism, or accurate extrapolation. The stored top-20 low-similarity list was a descriptive subset, not a validated set of PARG inhibitors.

## 5 Protein structures and receptor preparation

### 5.1 Reference system

The principal docking receptor was PARG chain A from PDB entry 6OAK, with the crystallographic M0M ligand defining the reference pocket. The laboratory record reports a resolution of 1.70 Å. The ligand reference contained 22 heavy atoms, and the docking centre was fixed at x=27.20, y=2.84, z=0.17 Å. Protein completeness and residues near the pocket were inspected before the original docking experiments. The existing prepared PDBQT files are the inputs to the current execution; the new run does not reconstruct their protonation or charge assignment.

The receptor-selection script retained chain A, removed M0M, and either removed all crystallographic waters or retained the five specified water residues. Other nonwater residues were accepted by the script. Inspection of the wet receptor PDB confirmed retained CME records as well as water records; the procedure should therefore not be described as removal of all nonprotein heteroatoms. The surviving preparation code does not establish a complete, explicit protocol for receptor protonation, pH-dependent states, partial charges, missing-residue modelling, or force-field minimization. Those details must be recovered from original preparation commands before a fully specified receptor-preparation claim can be made.

### 5.2 Operational water selection and its audit

In scripts/10_waters.py, each chain-A crystallographic water was represented by its first atom. A water was selected when it lay within 3.5 Å of at least one atom of M0M and within 3.5 Å of a protein N or O atom from a standard residue. Crucially, list(ligand) included deposited hydrogens. No donor–acceptor angle, ligand-atom chemical typing, occupancy cutoff, or B-factor cutoff was used. This operational rule selected HOH 1211, 1232, 1307, 1317, and 1402, which were then hard-coded into the wet-receptor preparation script.

The 16 October audit repeated the distance calculation both with all deposited ligand atoms and with ligand heavy atoms only. The all-atom calculation reproduced all five waters; the heavy-atom-only calculation retained only HOH 1317. Four waters were selected through proximity to hydrogens bonded to carbon. The wet condition is therefore described as retention of five geometrically selected crystallographic waters, without asserting that all five are conventional ligand–water–protein hydrogen-bond bridges. The previously written heavy-atom-only selection rule does not reproduce the implemented wet receptor.

| Water | Nearest ligand atom | All atom distance Å | Nearest heavy atom | Heavy atom distance Å | Protein N or O distance Å |
| --- | --- | --- | --- | --- | --- |
| 1211 | H121 | 3.006 | C12 | 3.800 | 2.754 |
| 1232 | H182 | 2.993 | C18 | 3.913 | 2.775 |
| 1307 | H182 | 2.770 | C18 | 3.581 | 2.890 |
| 1317 | H051 | 2.221 | N05 | 2.886 | 3.303 |
| 1402 | H141 | 3.181 | C14 | 3.987 | 3.241 |

All distances are in the untransformed 6OAK coordinate frame. The audit is saved as water_selection_audit.json in the current run directory. This finding does not change the identity of the receptor used in the historical experiments or in the ongoing reproduction.

### 5.3 Cross docking structures and coordinate alignment

The archived structure inventory explicitly requested 24 PDB entries: 4A0D, 4B1G, 4B1H, 4B1I, 4B1J, 5A7R, 5LHB, 6HH6, 6HMK, 6HML, 6HMM, 6HMN, 6O9X, 6O9Y, 6OA0, 6OA1, 6OA3, 6OAL, 7KFP, 7KG0, 7KG1, 7KG6, 7KG7, and 7KG8. The 6OAK reference was stored separately. This is a fixed historical inventory, not a claim to include every currently available PARG structure.

Chain-A Cα atoms were matched by residue-number and insertion-code fields, accepting blank or A alternative locations. A proper Kabsch rotation and translation were calculated from all matched Cα coordinates, without iterative outlier rejection. The same transform was applied to each extracted ligand to place its crystallographic coordinates in the 6OAK reference frame. Structures with fewer than 100 matched Cα atoms were skipped. The record reports 490–497 matches for the retained comparisons.

Heteroatom groups representing waters, common buffers, ions, modified residues, and other explicitly listed excluded residue names were removed from ligand consideration. Ligand groups containing 15–50 heavy atoms were retained for inventory. Phosphorus-containing groups were classified operationally as substrate-like and excluded from the inhibitor set. Candidate inhibitor centroids had to lie within 8 Å of the reference box centre. This generated the following 17 structure–ligand cases.

| Structure group | Retained structure and ligand identifiers | Cases |
| --- | --- | --- |
| 5LHB and 6HM series | 5LHB_6WJ; 6HMK_7JC; 6HML_73L; 6HMM_7JB; 6HMN_70J | 5 |
| 6O9 and 6OA series | 6O9X_M0S; 6O9Y_M0V; 6OA0_M1D; 6OA1_M0P; 6OA3_M0M; 6OAL_M0P | 6 |
| 7KF and 7KG series | 7KFP_WD7; 7KG0_WDA; 7KG1_WDD; 7KG6_WDG; 7KG7_WDJ; 7KG8_WDM | 6 |

These are 17 crystallographic cases rather than 17 independent chemical series. M0P occurs in two structures, and M0M is represented by a noncognate 6OA3 structure in cross docking. The separate 6OAK_M0M redocking experiment is not one of the 17 cases. Series relatedness and repeated chemical identities limit independence in the structural comparison.

### 5.4 Ligand preparation

Superposed crystallographic ligands were converted from PDB to PDBQT using Open Babel with -h [7]. The original conversion used the crystallographic coordinates and did not generate a new three-dimensional conformer from SMILES for these benchmark ligands. Heavy-atom counts in PDB and PDBQT were compared before a ligand was admitted to the target list. Preserving atom count is a preparation check, not independent verification of bond orders, tautomer choice, formal charge, or stereochemistry.

Library-control ligands were instead generated from SMILES with Open Babel --gen3d -h and written as PDBQT. An explicit pH, protonation-state enumeration procedure, and conformer-generation random seed were not supplied. These differences between crystal-derived benchmark inputs and generated library inputs are retained in the methodological record.

## 6 Docking experiments and pose evaluation

### 6.1 Vina settings

AutoDock Vina 1.2.5 was used with its default Vina scoring function, exhaustiveness 16, num_modes 9, and random seed 42 [8]. Receptors were rigid and ligand flexibility was defined by the PDBQT torsion tree. Up to nine output modes were requested; the number actually returned is recorded per job. The 24 Å cubic box was centred at (27.20, 2.84, 0.17) Å for the original dry/wet cross-docking experiments. The 16 Å control-screening box used the same centre. Other Vina settings, including the output energy window, were left at binary defaults. No flexible-side-chain docking, molecular dynamics, Boltz prediction, free-energy perturbation, or MM/PBSA calculation was performed in this workflow.

The historical native redocking placed M0M into the dry and wet 6OAK receptors. Cross docking placed each of the 17 crystallographic ligands into both receptors. The current fresh run repeats these operations with four Vina CPU threads, separate output paths, captured logs, and a 3,600-second timeout per docking call. A fixed seed does not make differently prepared ligand conformations or receptor files interchangeable.

### 6.2 Revision of the RMSD calculation

Sequential atom-order RMSD was initially found unsuitable because PDBQT conversion can reorder atoms. An intermediate method minimized the sum of squared interatomic distances under an element-identity constraint using linear_sum_assignment. This ignores molecular connectivity and supplies only a lower bound on a chemically constrained atom mapping. Results based on that lower bound were superseded for the main cross-docking analysis.

The revised calculation used the ideal ligand structure from the PDB Chemical Component Dictionary as a template. Hydrogen atoms were explicitly removed from reference and pose molecules. PDB molecules were loaded without sanitization, followed by property-cache updating and ring identification. The template was converted to a query with generic bonds after attempting kekulization and aromatic-flag clearing. Query matches were enumerated in both the crystallographic reference and the docked pose, with a limit of 300 matches per molecule in the historical scripts.

For each pair of full-template mappings, RMSD = sqrt[sum(||pose_i − reference_i||²) / N]. The minimum across the enumerated mapping pairs was retained. No additional rotation or translation of the ligand was fitted during RMSD evaluation: the receptor-based alignment had already established a common coordinate frame. The metric therefore evaluates pose placement rather than superposed ligand shape. It is best described as connectivity-constrained, symmetry-aware heavy-atom RMSD; the generic-bond query does not enforce bond order or stereochemistry and may permit mappings beyond strict chemical equivalence.

The historical top-pose script allowed a maximum-common-substructure fallback with a 20-second search timeout and a minimum of four atoms. The archived 34-row corrected table reports full-template matching for every row, so that fallback was not used in those results. The all-pose script used no MCS fallback. The new execution checks equality of template, reference, and pose heavy-atom counts, requires valid full-template matches, and probes with 301 matches to reject a mapping search that would exceed the historical cap. A preparation preflight passed for 53 stored top poses: two native redocks plus 17 ligands in each of the dry, wet, and B receptors. This preflight is not a fresh docking result.

### 6.3 Pose success sampling and ranking

A near-native pose was defined operationally as RMSD strictly below 2.0 Å. Top-1 success used the first Vina output mode; sampling success required at least one returned mode below the threshold. The crystallographic reference was used only for retrospective pose evaluation, not for score-based pose selection. All returned modes were evaluated so that generation of a near-native pose could be distinguished from its placement at the top of the ranked list.

The historical matched ranking analysis retained cases with a near-native pose in both dry and wet receptors and compared the rank of the minimum-RMSD pose. This quantity is not necessarily the rank of the earliest acceptable pose below 2.0 Å. These definitions must remain distinct in manuscript tables and figure labels. Threshold-sensitive cases, including the archived 6OA3_M0M wet value of 2.03 Å, are reported with their continuous RMSD rather than being assigned a mechanistic explanation solely from threshold crossing.

### 6.4 Second receptor and box diagnostics

A second receptor was prepared from 7KG1 chain A and aligned to 6OAK using matched Cα atoms. Its script selected waters within 3.5 Å of a WDD heavy atom and of protein N or O, reporting four retained waters. Unlike the 6OAK selection script, this ligand-distance calculation explicitly excluded hydrogens. The two receptor-water selection implementations therefore were not identical. Protein ATOM records and the selected waters were written into the aligned B-receptor PDB.

The 17 ligands were docked into receptor B with the original 24 Å box. Historical ensemble selection chose the receptor whose top-ranked pose had the more favourable Vina score, with ties assigned to receptor A. An oracle comparison selected the smaller retrospective RMSD between receptor-specific top poses; it was a diagnostic upper bound and was not an operational screening method. Historical ensemble RMSDs were calculated with the element-only lower bound and must not be pooled with corrected cross-docking RMSDs.

Diagnostic scripts examined receptor coordinate alignment, contacts between the reference ligand and protein, centroid displacement, atom-count changes during conversion, and waters around alternative poses. A tighter-box experiment used 16 Å cubes centred separately on each aligned reference-ligand centroid, rounded to two decimals. That condition uses knowledge of each crystallographic ligand position and is distinct from the fixed-centre library-screening protocol. The original record describes this 34-job experiment as interrupted after partial execution. The fresh runner queues all 34 jobs and retains their corrected all-pose RMSDs, but completion is determined from its logs rather than assumed here.

## 7 Structural analyses of water geometry and side chains

Binding-site residues were defined by at least one reference protein heavy atom lying within 6 Å of a native M0M heavy atom. After Cα alignment, matched side-chain heavy atoms were compared by residue identifier and atom name, excluding N, CA, C, O, and OXT and requiring matching residue names. Site side-chain RMSD was the square root of the mean squared displacement over those atom pairs. Maximum displacement was the largest individual matched-atom displacement in the site. Residue-specific displacement used the largest displacement within that residue. These coordinate comparisons did not explicitly correct symmetry-related exchanges of equivalent side-chain atom names.

The original diagnostic scripts defined their groups using best-pose RMSD from the obsolete lower-bound table, although some prose and figure labels described top-ranked docking success. These groupings and their aggregate displacement values must be regenerated from the chosen corrected outcome before reuse. The existing side-chain scatter plot is retained as historical output, not a corrected figure for the new run.

Tyr795 and Phe902 χ1 angles were calculated from the N–CA–CB–CG coordinates using an atan2 dihedral expression. Side-chain B factors were averaged over the included atom records as a structural descriptor. Rotamer assignments described two observed coordinate groups. Angular separation must be interpreted on a circle; the shortest separation between the reported Tyr795 means of +66.3° and −168.5° is approximately 125.2°, not 234.8°. Static structures and B factors do not establish exchange kinetics or prove that an individual residue dynamically samples both states.

For water-contact analysis, the minimum distance between an aligned crystallographic ligand heavy atom and each retained 6OAK water oxygen was measured. These cross-structure distances assess compatibility with the fixed reference-water positions. They do not describe a water actually deposited with every ligand. A short contact may indicate geometric overlap; a distance alone, particularly to carbon or sulfur, does not establish a hydrogen bond. The original binary water-bridging versus water-displacing classification was exploratory. Independent alignment of 7KFP placed WD7 approximately 0.94 Å from HOH 1402, demonstrating a conflict with its original blanket series assignment. Class-sensitive analyses require explicit per-water distances and alternative-assignment sensitivity reporting.

## 8 Random library control and failed jobs

The 100 highest-ranked library entries were compared with 100 randomly selected entries from the same library after excluding the top-100 identifiers. Python random.seed(42), followed by random.sample, defined the historical control sample. Jobs were interleaved in top/random pairs. Both groups used the wet 6OAK receptor, the fixed-centre 16 Å box, exhaustiveness 16, seed 42, and up to nine modes. The outcome was the top-ranked Vina score in kcal/mol; lower values were considered more favourable within this docking protocol.

The original scripts reused ligand and docking files whenever they already existed and did not associate cache reuse with hashes of receptor, ligand, and settings. They also suppressed subprocess output and saved successful docking records only. The historical CSV contained 92 top-ranked and 100 random outcomes, with eight missing top-ranked entries. Those entries were documented as closely related stereoisomeric forms. The original logs do not establish whether each missing outcome arose in preparation, search, or another execution step; failure is not evidence of weak binding.

The fresh run starts with a separate ligand/output directory, captures preparation and Vina logs, checks return codes, and imposes a 600-second preparation timeout and a 3,600-second docking timeout. All 200 attempted control entries are written to control_all.csv, with preparation return codes and statuses of ok, preparation_failed, or docking_failed. Timeout return code 124 is recorded. Failure remains a missing scientific outcome even when its computational cause is known.

## 9 Statistical analyses

Historical continuous-score comparisons used the Mann–Whitney U test. The original control script used alternative='less' for the top-ranked group; the independent audit additionally calculated the two-sided value. Enrichment below −9 kcal/mol was evaluated with Fisher's exact test, originally using alternative='greater' and subsequently also reporting the two-sided test. Means, medians, extrema, threshold counts, and attempted-versus-successful denominators were retained. Docking-score enrichment is not enrichment in experimentally confirmed PARG inhibitors.

Sensitivity analyses replaced the eight missing top-group scores with −6.42 kcal/mol, the least favourable score observed across the two groups, or with −10.21 kcal/mol, the most favourable observed top-group score. The second value was not the best score in the full control experiment, which included a random-group score of −10.60 kcal/mol. These are specified imputation scenarios, not confidence limits or a model of the missing-data mechanism. The fresh runner's automated report currently supplies a two-sided complete-case Mann–Whitney comparison and failure counts; it does not automatically reproduce every historical imputation or enrichment analysis.

Paired dry/wet binary pose outcomes were assessed using the exact McNemar test, equivalent to a two-sided binomial test with probability 0.5 on discordant pairs. Paired continuous RMSDs and the restricted matched ranking outcomes were assessed historically with Wilcoxon signed-rank tests. The class comparison used the difference wet RMSD minus dry RMSD and a between-class Mann–Whitney test. Sampling failure to reach statistical significance was not interpreted as equivalence. At six matched cases, an all-same-direction exact signed-rank comparison can attain a two-sided p value of 0.03125, illustrating the coarse inferential resolution of the small sample.

Analyses were exploratory. No multiple-comparison correction, preregistered primary hypothesis, scaffold-cluster resampling, or experimental validation was implemented. Ligands within medicinal chemistry series and repeated chemical identities are not fully independent. Nominal p values therefore accompany their denominators and effect descriptions rather than supporting broad mechanistic or therapeutic claims. All new automated summary p values are calculated from the new outputs rather than copied from historical figures.

## 10 Reproducibility and data management

Calculations were performed in Ubuntu 24.04 under WSL2 using the parg conda environment. The inspected runtime reported eight logical CPUs, approximately 7.7 GiB of WSL memory, and no available CUDA GPU. The PyTorch package contained CUDA support, but that package property did not imply GPU use. The exact imported environment was recorded separately from environment.yml and requirements.txt.

| Software | Recorded runtime version | Main role |
| --- | --- | --- |
| Python | 3.11.16 | Pipeline and analysis |
| RDKit | 2026.3.6 | Structures fingerprints descriptors |
| pandas | 3.0.6 | Tabular data |
| NumPy | 2.4.6 | Numerical arrays |
| SciPy | 1.17.1 | Statistics and assignment |
| scikit learn | 1.9.1 | Evaluation metrics |
| LightGBM | 4.7.0 | Gradient boosting |
| Chemprop | 2.3.1 | Graph neural network |
| PyTorch | 2.14.1 with cu130 build | Neural network backend |
| Open Babel | 3.2.1 | Structure conversion and generation |
| AutoDock Vina | 1.2.5 | Docking |

The fresh run copies the raw activity table, library snapshot, CCD templates, prepared structures, analysis scripts, fixed split, and saved Chemprop configurations to a dated directory. Input files in the data, structures, scripts, and configs directories receive SHA-256 hashes in input_manifest.json. The Vina executable is copied into the run directory; the initial manifest generator does not include that root-level binary, so its checksum should be recorded separately when packaging the run. Original project outputs are preserved.

The runner records stage names, return codes, elapsed times, and the current state in status.json. Each stage has its own log. Benchmark predictions, metrics, final feature names, screening failures, docking commands, and control failures are retained. The scoring pipeline was checked using 24 valid molecules and one deliberately invalid SMILES; all valid predictions matched direct model predictions and the invalid entry was recorded. These execution checks validate the implementation paths tested, not model applicability or docking chemistry.

The historical figure scripts mixed direct CSV inputs with transcribed numerical arrays and fallback values. The fresh runner generates its model benchmark and paired cross-docking plot directly from current CSV outputs. It does not yet regenerate every historical rotamer, side-chain, pose-panel, sensitivity, or ensemble figure. Those outputs retain their individual provenance and correction status.

## Appendix A Execution status and verified numerical anchors

The status snapshot below is populated from the run directory when this document is built. It is a timestamped record, not a continuously updating dashboard. A completed historical analysis and a newly completed reproduction are separate pieces of evidence.

Status captured 2026-10-16T13:34:03+04:00. Run state: running. Current stage: 06_screen. Completed stages: 01_curate, 02_split, 03_ml, 04_chemprop_seed0, 04_chemprop_seed1, 04_chemprop_seed2, 05_chemprop_evaluation. At this snapshot, screening had processed 110,592 entries and scored 110,592 successfully.

The fresh model benchmark reproduced scaffold-test RMSE 0.509188 and R² 0.560629 for the five-neighbour model, and RMSE 0.445215 and R² 0.664097 for LightGBM. The training-mean predictor yielded RMSE 0.800719 and R² −0.086514. The negative mean-baseline R² replaces the hard-coded zero in earlier summaries. On the random split, the corresponding RMSE values were 0.411186, 0.338293, and 0.738076. These are measured outputs of the completed new benchmark stages.

The historical corrected cross-docking tables contained 34 top-pose cases and 306 poses overall. They reported top-1 success in 0/17 dry and 6/17 wet cases, with a near-native pose among returned modes in 6/17 dry and 7/17 wet cases. These historical figures must not be presented as outcomes of the fresh docking stage until that stage completes. Archived native-redocking, ensemble, and some side-chain summaries still used the older lower-bound metric and are not substituted for corrected new outputs.

## Appendix B File and script provenance

All relative paths in this appendix are relative to /home/aliwhoami/parg unless prefixed by the dated run directory.

| Analysis | Primary implementation or record | Principal output |
| --- | --- | --- |
| Curation and splits | scripts/02_curate.py and 03_split.py | data/parg_curated.csv and parg_split.csv |
| Historical baseline and boosting | scripts/04_baseline.py and 05_lightgbm.py | Terminal logs and model comparison records |
| Chemprop | models/chemprop_seed0 through seed2/config.toml; data/fixed_split.json | best.pt and test_predictions.csv |
| Final model and screen | scripts/29_train_final.py; 30_get_zinc.py; 32_screen.py; 33_novelty.py | lgbm_final.pkl and screen_ranked files |
| Water and receptor selection | scripts/10_waters.py and 11_prepare_receptor.py | Receptor PDB and prepared PDBQT inputs |
| Structure inventory and alignment | scripts/15_inventory.py through 17_prep_ligands.py | structures/xdock and xdock_targets.txt |
| Historical docking | scripts/18_crossdock.py and 24_ensemble.py | PDBQT poses and old metric tables |
| Revised RMSD | scripts/38_true_rmsd.py and 40_allposes.py | crossdock_truermsd.csv and crossdock_allposes.csv |
| Side chains and rotamers | scripts/19_diagnose.py; 21_which_residue.py; 22_rotamer.py | Historical terminal summaries |
| Alternative receptor diagnostics | scripts/23_build_stateB.py through 28_box16.py | B receptor and diagnostic poses |
| Random controls and sensitivity | scripts/34_control.py and 35_sensitivity.py | control.csv and historical statistics |
| Fresh execution | runs/20261016_132012/run_parg_project.py | status.json; logs; new models and results |
| Water selection audit | runs/20261016_132012/water_selection_audit.json | All atom versus heavy atom distances |

## Appendix C Reporting boundaries

The activity model pooled IC50, EC50, and Kd records; it should not be called an assay-homogeneous biochemical affinity model. The five-water receptor derives from an all-ligand-atom selection rule; it should not be attributed to the unreproduced heavy-atom rule. The B-receptor water-selection implementation differs from the A-receptor implementation. Prepared structures are reused in the current run, with their remaining chemical-preparation uncertainties preserved rather than silently filled in.

The current RMSD method constrains connectivity but not every aspect of chemical equivalence. A full atom count and successful query match do not prove correct protonation or stereochemistry. Historical comparisons based on minimum-RMSD pose rank, earliest acceptable pose rank, top-1 success, and sampling success must be labelled separately. The small structural dataset, chemical-series dependence, missing control outcomes, and exploratory testing limit inference.

No completed evidence was identified for molecular dynamics, Boltz screening, experimental assays, an independently validated uncertainty model, prospective supplier availability, or a public reproducibility deposition. These activities must not be written in the past tense as completed methods. The present document records the implemented work and the explicit execution changes, including diagnostic inconsistencies that affect reproducibility.

## References

1. Zdrazil B, et al. The ChEMBL Database in 2023 a drug discovery platform spanning multiple bioactivity data types and time periods. Nucleic Acids Research. 2024;52(D1):D1180–D1192. https://doi.org/10.1093/nar/gkad1004

2. Rogers D, Hahn M. Extended connectivity fingerprints. Journal of Chemical Information and Modeling. 2010;50(5):742–754. https://doi.org/10.1021/ci100050t

3. Ke G, et al. LightGBM A Highly Efficient Gradient Boosting Decision Tree. Advances in Neural Information Processing Systems. 2017;30. https://proceedings.neurips.cc/paper_files/paper/2017/hash/6449f44a102fde848669bdd9eb6b76fa-Abstract.html

4. Heid E, et al. Chemprop A Machine Learning Package for Chemical Property Prediction. Journal of Chemical Information and Modeling. Published online December 26, 2023. https://doi.org/10.1021/acs.jcim.3c01250

5. Chemprop v2 An Efficient Modular Machine Learning Package for Chemical Property Prediction. ChemRxiv preprint record, 2025. https://doi.org/10.26434/chemrxiv-2025-4p1nr

6. Irwin JJ, et al. ZINC20 A Free Ultralarge Scale Chemical Database for Ligand Discovery. Journal of Chemical Information and Modeling. 2020;60(12):6065–6073. https://doi.org/10.1021/acs.jcim.0c00675

7. O'Boyle NM, et al. Open Babel An open chemical toolbox. Journal of Cheminformatics. 2011;3:33. https://doi.org/10.1186/1758-2946-3-33

8. Eberhardt J, Santos Martins D, Tillack AF, Forli S. AutoDock Vina 1.2.0 New Docking Methods Expanded Force Field and Python Bindings. Journal of Chemical Information and Modeling. 2021;61(8):3891–3898. https://doi.org/10.1021/acs.jcim.1c00203

Software citations identify the underlying tools; exact installed versions and implementation details are given in the methods and runtime record. The Chemprop v2 source above is identified as a preprint rather than assumed to be a peer-reviewed final publication.
