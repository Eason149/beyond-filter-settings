# Beyond filter settings: How viewers judge edited faces

Analysis and plotting code for the study of editing operations, Geometry–Appearance–Identity (GAI) descriptors, beauty and perceived unnaturalness ratings, and viewing-period EEG. The repository distributes source code and [aggregate research results](data/aggregate/). Human-participant records, facial images, and model weights are excluded.

## Data and code availability

The analysis and visualization code is available at [commit `2cc32c5`](https://github.com/Eason149/beyond-filter-settings/tree/2cc32c5417d79a931ebf28fc7449a853fd12831c).

The manuscript and supplementary material report aggregate estimates. Six downloadable CSV tables report behavioral and EEG contrasts, appearance-model estimates and comparisons, bootstrap intervals, and held-out prediction metrics. See the [aggregate data directory at commit `52dfc9e`](https://github.com/Eason149/beyond-filter-settings/tree/52dfc9eb6add64e82b5cb5c9e83dc45565691fef/data/aggregate) for the tables and their descriptions.

Participant-level EEG recordings, trial-level behavioral records, facial stimuli, and landmark coordinates are not distributed in this repository. Access to those materials is subject to the applicable participant consent and image permissions. See [DATA_ACCESS.md](DATA_ACCESS.md) for the input layout and access notes. For data access requests, please contact the author, Yi Cao, at [cy1840156@gmail.com](mailto:cy1840156@gmail.com) or [cyi06@vt.edu](mailto:cyi06@vt.edu).

## Code map

| Manuscript component | Entry points |
|---|---|
| Source measurements and EEG preprocessing | `scripts/metrics/`, `scripts/preprocessing/` |
| Operations and rating/EEG contrasts | `scripts/analysis/run_stage28_ipm_route_search.py` |
| Realized visual changes and ratings | `run_stage26_frozen_metric_models.py` and `run_stage27_a2_validated_route.py` under `scripts/analysis/` |
| GAI–EEG and EEG–evaluation models | `scripts/analysis/run_stage29_cascade_validation.py` |
| Held-out prediction, component comparisons, identity omission, crossed mixed models and EEG sensitivity | `ipm_revision_round1/scripts/04_prediction.py`, `05_mixed_models.R`, `06_eeg_sensitivity.py` |
| Exploratory whole-scalp cluster analysis underlying the channel–time display | `scripts/analysis/run_paper_grade_extension.py` |
| Statistical figure layouts and single-column variants | `plotting/source/` |

## Running with authorized inputs

Use the repository-level environment files for the historical Python analyses and `ipm_revision_round1/requirements.txt` for the later prediction and sensitivity models. The crossed models require R and **lme4**. Set `IPM_DATA_ROOT` to the authorized study data tree and `IPM_STIMULUS_ROOT` for the protected face images. See [DATA_ACCESS.md](DATA_ACCESS.md). Baseline reproduction scripts also read previously derived reference estimates. Existing archive scripts retain their historical stage identifiers.

From the repository root, authorized users with the private inputs can run `python ipm_revision_round1/scripts/04_prediction.py`, `Rscript ipm_revision_round1/scripts/05_mixed_models.R ipm_revision_round1`, and `python ipm_revision_round1/scripts/06_eeg_sensitivity.py`. The R model expects the prepared `ipm_revision_round1/tables/behavior_input.csv`; the Python loaders document the required columns and input layout. A clean checkout alone cannot reproduce the numerical results.

For the R figures, supply the private figure-input CSVs named in `draw_figures.R` under `$IPM_FIGURE_ROOT/data/` (default: `plotting/data/`). The figures use ggplot2, patchwork, jsonlite, svglite, ragg, png and qpdf; the column layouts additionally use pdftools for PDF previews. Run the complete `draw_figures.R` before the column-specific scripts. The Python interpreter for layout checks is configurable with `FIGURE_QA_PYTHON`. No target images or participant files are published alongside the code.

## License

The adapted panel-alignment utilities are accompanied by their upstream license in `plotting/source/LICENSE-nature-skills.txt`. No additional reuse license for the authors’ code or protected data is granted by this release.


