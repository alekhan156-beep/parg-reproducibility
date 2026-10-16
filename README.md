# PARG computational screening: reproducibility snapshot

Prepared for research review on 2026-10-16. This is an interim snapshot, not a final publication release.

## Contents and provenance

- `scripts/`, `structures/`, `models/`, and small `data/` files preserve the supplied project. Root models are historical models.
- `snapshot_20261016/` contains freshly trained models, completed screening and benchmark tables, configurations, input hashes and captured run status from run `20261016_132012`.
- `historical/` preserves earlier results, figures, desktop scripts and lab records. These can disagree with the fresh run and must not be combined as one experiment.
- `docs/PARG_Methodology.pdf` and its Markdown source document methods and limitations as of their stated snapshot time.
- The exact executed runner is in `snapshot_20261016/run_parg_project.py`. The root runner changes only source-directory and Python-environment path resolution for portability. This portable version has been syntax-checked, not run end to end.

## Reproduction in Ubuntu / WSL

1. Install Conda, Git and the system dependencies recorded by the project. Create the environment with `conda env create -f environment.yml`, then activate `parg`. The environment export is platform-specific and has not been rebuilt on a clean machine.
2. Extract the companion `parg-large-inputs.tar.gz` into this repository. It supplies the frozen ZINC library and the exact Linux Vina executable. Run `chmod +x vina_1.2.5_linux_x86_64`.
3. Run `python run_parg_project.py` from this directory. This launches the full workflow in a new `runs/` directory and may take hours. Read `runs/LATEST`, the run's `status.json`, and stage logs to check completion.
4. Compare generated tables with `snapshot_20261016/results/`. Full historical screening scores are in the companion archive; newly generated outputs remain in the new run directory.

Historical scripts can contain original absolute paths. Use the root runner for the documented rerun. Neither historical scripts nor this README's commands have been independently validated on another computer.

## Scientific status and limitations

At packaging, curation, ML training, three Chemprop runs, full-library screening and novelty analysis had completed. Docking/control analysis was still running; in-flight control files are excluded. Final statistics and figures require a later release.

The final screening predictor is a LightGBM regressor trained on all 870 curated compounds. The dataset combines IC50, EC50 and one Kd observation before compound aggregation; these endpoints are not biologically interchangeable. Random Forest has not been benchmarked in this workflow. Predictions and docking scores do not establish experimental activity.

Prepared receptor chemistry is reused. The five-water selection used all deposited ligand atoms, including hydrogens; it does not implement the stated heavy-atom-only rule. Current pose RMSD matching uses connectivity and does not fully enforce bond order or stereochemistry. The small cross-docking set contains related ligands. Historical figures may use different RMSD calculations and must be checked before publication. See the methodology and water audit for details.

## Files, access and licensing

`MANIFEST.sha256.json` records hashes of packaged files. Companion archives have separate hashes in `release_assets/release_manifest.json`. No credentials or environment secrets are intentionally included; public redistribution still requires a provenance/license review. No new software license is granted by this snapshot. Third-party data and executable rights remain with their respective owners. The owner explicitly authorized public upload of this interim snapshot on 2026-10-16.

## Reconstruct the companion archive

The frozen library, complete screening scores, and exact Linux Vina executable are included in `release_assets/` as numbered archive parts to keep each Git file below the per-file size limit. These are data files, not Git LFS pointers.

From the repository root on Ubuntu/WSL:

```bash
cat release_assets/parg-large-inputs.tar.gz.part* > parg-large-inputs.tar.gz
sha256sum parg-large-inputs.tar.gz
tar -xzf parg-large-inputs.tar.gz
chmod +x vina_1.2.5_linux_x86_64
```

Expected archive SHA-256: `899a9724069e7d16d5e426cdbb87f11806bd0344f200fc3f7d2c008fcd4eadef`.

The repository contains binary model files and compressed data; cloning therefore requires several hundred megabytes. The original package ZIP listed in the release manifest is a local convenience copy, not a required input.
