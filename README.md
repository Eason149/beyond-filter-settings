# IPM Facial Editing Analysis Code

This public repository contains **analysis code only** for a behavioral and EEG study of algorithmic facial editing. It is separated from a private archival repository that retains the protected research inputs, frozen derivatives, result tables, figures, audit reports, and manuscript materials.

## What is public

- Python scripts for EEG preprocessing and ERP analysis;
- behavioral, EEG--behavior, sensitivity, and specification-audit workflows;
- visual-information metric code for geometric change (`G`) and appearance change (`A`/`A2`);
- MVPA and representational-similarity analysis workflows;
- Python/Conda dependency specifications and path templates.

## What remains private

- raw or processed participant EEG and behavioral data;
- original or edited face images and facial-landmark coordinate files;
- participant-level and aggregate numerical results;
- generated tables, figures, logs, model artifacts, and route-decision reports;
- the manuscript, supplements, internal audits, and working notes.

This split is deliberate. The public repository supports inspection of the computational specifications without publishing human-participant data, facial materials, unpublished numerical findings, or manuscript content. The absence of data and outputs means that a fresh clone is **not independently executable** and should not be interpreted as a data release or as evidence for any empirical claim. See [DATA_ACCESS.md](DATA_ACCESS.md).

## Repository structure

```text
scripts/
  analysis/       behavioral, EEG, joint-model, sensitivity, and route analyses
  metrics/        G, A/A2, landmark, and image-input validation
  preprocessing/  Curry-to-epoch preprocessing and ERP/MVPA analyses
  rsa/            time-resolved and component-level RSA workflows
data/README.md    placeholder and disclosure boundary; no data are included
environment.yml  Conda environment specification
requirements.txt pip dependency specification
```

## Environment

The reference environment uses Python 3.11. With Conda:

```bash
conda env create -f environment.yml
conda activate ipm-face-editing-eeg
```

Alternatively, create an isolated Python 3.11 environment and run:

```bash
python -m pip install -r requirements.txt
```

MediaPipe is fixed at `0.10.21` for the 468-landmark reproduction workflow. Exact random seeds, permutation counts, ROIs, time ranges, missing-trial handling, and correction families are defined in the relevant scripts.

## Private input configuration

Do not copy protected inputs into this Git repository. Point the code to authorized local storage instead:

```powershell
$env:IPM_DATA_ROOT = 'X:\authorized\project_data'
$env:IPM_STIMULUS_ROOT = 'X:\authorized\core_images'
```

On macOS/Linux:

```bash
export IPM_DATA_ROOT=/authorized/project_data
export IPM_STIMULUS_ROOT=/authorized/core_images
```

Several later-stage scripts consume outputs produced by earlier stages. Their required schemas and paths are visible in the source. A reproducibility run should preserve the declared sample rules and analysis specifications; it must not select samples, metrics, identities, or time windows according to statistical significance.

## Main workflow map

1. Inspect and preprocess authorized Curry recordings with `scripts/preprocessing/preprocess_raw_curry.py`.
2. Run primary and secondary ERP workflows in `scripts/preprocessing/`.
3. Audit and compute visual metrics with `scripts/metrics/`.
4. Reconcile the EEG sample and run frozen-metric models with the Stage 2.6 scripts in `scripts/analysis/`.
5. Run the validated A2 and stimulus-bounded sensitivity workflows with the Stage 2.7--2.8 scripts.
6. Treat MVPA, RSA, whole-scalp, and EEG--behavior analyses according to their multiplicity families and stated inferential boundaries.
7. Keep manuscript-specific figures and result narratives in the private archive.

The code contains exploratory and multiplicity-controlled analyses. It does not by itself establish confirmatory mechanisms, causal pathways, universal face-processing mechanisms, or generalization to new identities.

## Citation and license

Citation metadata are provided in `CITATION.cff`. No open-source license has been assigned yet; public visibility does not waive copyright or grant permission to reuse the code beyond rights provided by law. A license can be added after the authors select one.

## Private archive

A separate private repository is maintained by the project owner as the controlled archival record. It contains the non-public analysis outputs and writing materials and is not converted to public visibility by this code release.
