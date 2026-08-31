# Data access and disclosure boundary

This is a **code-only public repository**. It contains no participant data, facial stimuli, participant-level derivatives, aggregate result tables, generated figures, manuscript files, or result narratives.

The protected inputs include Curry EEG recordings, processed EEG epochs, E-Prime logs, trial-level behavioral observations, original and edited facial images, 468-point facial landmarks, image/model hashes that could reveal protected materials, and trained model binaries. These materials are retained outside this public repository under the applicable consent, institutional, licensing, and storage restrictions.

The expected private working layout is documented only to make path dependencies explicit:

```text
${IPM_DATA_ROOT}/
  EEGDATA/EEGDATA/curry8/Daya-zqy26/sXX/*.cdt
  EEGDATA/EEGDATA/eprime/
  CHB_multimodal_facial_editing/
  paper_extension_final/reanalysis_30/
  erp_dynamics_final/
  ipm_visual_information_gate/
  ipm_stage_2_6/
  ipm_stage_2_7_a2_validated_route/
  mvpa_final/
  RSA_multimodal_geometry/
```

`IPM_STIMULUS_ROOT` should point to an authorized local directory containing the 68 core images. Filenames must match the private condition mapping used by the analysis.

Cloning this repository alone is therefore insufficient to reproduce numerical findings. Authorized researchers must obtain the protected inputs through the data custodian and apply the same directory mapping, versioned environment, sample rules, random seeds, and multiplicity corrections encoded in the scripts.

No protected data should be added to a fork, issue, pull request, Actions artifact, or public release. The repository `.gitignore` provides an additional safeguard but is not a substitute for checking every staged file before publication.
