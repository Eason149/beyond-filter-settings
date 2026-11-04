#!/usr/bin/env python
"""Paper-grade extension analyses for the edited-face EEG project.

This script complements, but does not replace, the confirmed ERP dynamics
analysis. It adds:
1. whole-scalp spatiotemporal cluster permutation tests;
2. formal stage-by-edit-class and rating-by-edit-class interactions;
3. temporal-centroid and leave-one-subject-out stability analyses;
4. ROI and temporal-smoothing robustness checks;
5. integrated paper-ready figures and machine-readable summaries.
"""

from __future__ import annotations

import json
import logging
import math
import os
import subprocess
import sys
from collections import deque
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np
import pandas as pd
from scipy import stats

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from run_eeg_behavior_integration import clean_behavior, subject_behavior_betas
from run_erp_dynamics_persistence import CHANNEL_POS, STAGES, design_matrix
from run_time_resolved_factor_decoding import (
    FACTORS,
    N170_LPP_CHANNELS,
    POSTERIOR_CHANNELS,
    STRUCTURAL,
    SURFACE,
    Context,
    extract_subject_matrix,
    fdr_bh,
    find_eeg_files,
    load_subject,
    natural_subject_key,
    parse_eprime_condition_map,
    parse_eprime_records,
)


CONTRASTS = [
    "FSlim",
    "Eye",
    "Mouth",
    "Skin",
    "Surface",
    "Structural",
    "Surface_minus_Structural",
]
PLOT_CONTRASTS = ["Surface", "Structural", "Surface_minus_Structural", "Skin", "Mouth"]
COLORS = {
    "FSlim": "#167d6d",
    "Eye": "#d66b27",
    "Mouth": "#6667ab",
    "Skin": "#b83280",
    "Surface": "#bb6b12",
    "Structural": "#087f8c",
    "Surface_minus_Structural": "#374151",
}
STAGE_ORDER = list(STAGES)


def setup_dirs(root: Path) -> dict[str, Path]:
    dirs = {
        "root": root,
        "tables": root / "tables",
        "figures": root / "figures" / "paper_ready",
        "supplement": root / "figures" / "supplement",
        "summaries": root / "summaries",
        "intermediate": root / "intermediate",
        "logs": root / "logs",
    }
    for path in dirs.values():
        path.mkdir(parents=True, exist_ok=True)
    return dirs


def setup_logging(path: Path) -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler(path, encoding="utf-8"), logging.StreamHandler(sys.stdout)],
        force=True,
    )


def save_csv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, index=False, encoding="utf-8-sig")
    logging.info("Wrote %s shape=%s", path, df.shape)


def save_figure(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def json_default(value: Any) -> Any:
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return None if np.isnan(value) else float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def holm(values: pd.Series | np.ndarray) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    out = np.full(len(p), np.nan)
    valid = np.where(np.isfinite(p))[0]
    if not len(valid):
        return out
    order = valid[np.argsort(p[valid])]
    adjusted = np.maximum.accumulate((len(order) - np.arange(len(order))) * p[order])
    out[order] = np.minimum(adjusted, 1.0)
    return out


def t_values(data: np.ndarray) -> np.ndarray:
    n = np.sum(np.isfinite(data), axis=0)
    mean = np.nanmean(data, axis=0)
    sd = np.nanstd(data, axis=0, ddof=1)
    denom = sd / np.sqrt(np.maximum(n, 1))
    return np.divide(mean, denom, out=np.zeros_like(mean), where=denom > 0)


def build_adjacency(channels: list[str]) -> list[set[int]]:
    positions = np.array([CHANNEL_POS[ch.upper()] for ch in channels], dtype=float)
    distances = np.sqrt(((positions[:, None, :] - positions[None, :, :]) ** 2).sum(axis=2))
    adjacency = [set() for _ in channels]
    for i in range(len(channels)):
        nearest = np.argsort(distances[i])[1:5]
        for j in nearest:
            if distances[i, j] <= 0.58:
                adjacency[i].add(int(j))
                adjacency[int(j)].add(i)
    return adjacency


def connected_clusters(mask: np.ndarray, adjacency: list[set[int]]) -> list[np.ndarray]:
    n_times, n_channels = mask.shape
    seen = np.zeros_like(mask, dtype=bool)
    clusters: list[np.ndarray] = []
    for time_idx, channel_idx in np.argwhere(mask):
        if seen[time_idx, channel_idx]:
            continue
        queue = deque([(int(time_idx), int(channel_idx))])
        seen[time_idx, channel_idx] = True
        cells = []
        while queue:
            t_idx, ch_idx = queue.popleft()
            cells.append((t_idx, ch_idx))
            for next_t in (t_idx - 1, t_idx + 1):
                if 0 <= next_t < n_times and mask[next_t, ch_idx] and not seen[next_t, ch_idx]:
                    seen[next_t, ch_idx] = True
                    queue.append((next_t, ch_idx))
            for next_ch in adjacency[ch_idx]:
                if mask[t_idx, next_ch] and not seen[t_idx, next_ch]:
                    seen[t_idx, next_ch] = True
                    queue.append((t_idx, next_ch))
        clusters.append(np.asarray(cells, dtype=int))
    return clusters


def find_signed_spatiotemporal_clusters(
    tvals: np.ndarray, threshold: float, adjacency: list[set[int]]
) -> list[tuple[str, np.ndarray, float]]:
    found = []
    for direction, mask in (("positive", tvals >= threshold), ("negative", tvals <= -threshold)):
        for cells in connected_clusters(mask, adjacency):
            mass = float(np.sum(np.abs(tvals[cells[:, 0], cells[:, 1]])))
            found.append((direction, cells, mass))
    return found


def spatiotemporal_cluster_test(
    data: np.ndarray,
    times: np.ndarray,
    channels: list[str],
    adjacency: list[set[int]],
    n_perm: int,
    seed: int,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    t_obs = t_values(data)
    threshold = float(stats.t.ppf(0.975, data.shape[0] - 1))
    observed = find_signed_spatiotemporal_clusters(t_obs, threshold, adjacency)
    null_max = np.zeros(n_perm, dtype=float)
    for pidx in range(n_perm):
        signs = rng.choice([-1.0, 1.0], size=(data.shape[0], 1, 1))
        perm_t = t_values(data * signs)
        clusters = find_signed_spatiotemporal_clusters(perm_t, threshold, adjacency)
        null_max[pidx] = max((item[2] for item in clusters), default=0.0)

    rows = []
    cluster_labels = np.zeros(t_obs.shape, dtype=int)
    for cluster_id, (direction, cells, mass) in enumerate(observed, start=1):
        p_value = float((np.sum(null_max >= mass) + 1) / (n_perm + 1))
        cell_t = np.abs(t_obs[cells[:, 0], cells[:, 1]])
        peak_cell = cells[int(np.argmax(cell_t))]
        used_channels = sorted({channels[idx] for idx in cells[:, 1]})
        rows.append(
            {
                "cluster_id": cluster_id,
                "direction": direction,
                "cluster_start_ms": float(times[cells[:, 0]].min()),
                "cluster_end_ms": float(times[cells[:, 0]].max()),
                "duration_span_ms": float(times[cells[:, 0]].max() - times[cells[:, 0]].min() + np.median(np.diff(times))),
                "n_cells": int(len(cells)),
                "n_channels": int(len(used_channels)),
                "channels": ";".join(used_channels),
                "cluster_mass_abs_t": mass,
                "cluster_p": p_value,
                "peak_time_ms": float(times[peak_cell[0]]),
                "peak_channel": channels[int(peak_cell[1])],
                "peak_t": float(t_obs[peak_cell[0], peak_cell[1]]),
            }
        )
        if p_value < 0.05:
            cluster_labels[cells[:, 0], cells[:, 1]] = cluster_id
    return pd.DataFrame(rows), t_obs, cluster_labels


def moving_average(data: np.ndarray, width: int) -> np.ndarray:
    if width <= 1:
        return data
    kernel = np.ones(width, dtype=float) / width
    return np.apply_along_axis(lambda x: np.convolve(x, kernel, mode="same"), 1, data)


def load_channel_betas(
    project_root: Path,
    times: np.ndarray,
    min_trials_per_cell: int = 8,
) -> tuple[np.ndarray, list[str], list[str], pd.DataFrame]:
    dummy_dirs = setup_dirs(project_root / "paper_extension_final" / "_loader")
    ctx = Context(
        project_root=project_root,
        outdir=dummy_dirs["root"],
        seed=20260731,
        n_permutations=0,
        quick_test=False,
        dirs={
            "root": dummy_dirs["root"],
            "tables": dummy_dirs["tables"],
            "figures": dummy_dirs["root"],
            "paper_ready": dummy_dirs["root"],
            "logs": dummy_dirs["logs"],
            "summaries": dummy_dirs["summaries"],
            "intermediate": dummy_dirs["intermediate"],
            "tgm_matrices": dummy_dirs["intermediate"],
            "models": dummy_dirs["intermediate"],
        },
    )
    condition_map = parse_eprime_condition_map(project_root, ctx)
    subjects = [load_subject(path, ctx, condition_map) for path in find_eeg_files(project_root)]
    subjects = [subject for subject in subjects if subject is not None]
    common = set(subjects[0].channels)
    for subject in subjects[1:]:
        common &= set(subject.channels)
    channels = [ch for ch in subjects[0].channels if ch in common and ch.upper() in CHANNEL_POS]
    logging.info("Whole-scalp analysis uses %d common positioned channels: %s", len(channels), channels)

    subject_arrays = []
    used_subjects = []
    qc_rows = []
    half_width = float(np.median(np.diff(times)))
    for subject in subjects:
        meta = subject.metadata[
            (~subject.metadata["is_control"])
            & (~subject.metadata["attention_check"])
            & subject.metadata["raw_cond_id"].isin(range(2, 18))
        ].copy()
        counts = meta.groupby(FACTORS).size()
        qc = {
            "subj": subject.subj,
            "n_trials": int(len(meta)),
            "n_factor_cells": int(len(counts)),
            "min_trials_per_cell": int(counts.min()) if len(counts) else 0,
            "included": False,
            "reason": "",
        }
        if len(counts) < 16 or counts.min() < min_trials_per_cell:
            qc["reason"] = "incomplete factorial cells or insufficient cell count"
            qc_rows.append(qc)
            continue
        lookup = {ch: idx for idx, ch in enumerate(subject.channels)}
        channel_idx = np.array([lookup[ch] for ch in channels], dtype=int)
        subject.channel_indices = channel_idx
        trial_data = extract_subject_matrix(subject, meta, subject.times)
        windowed = []
        for center in times:
            mask = (subject.times >= center - half_width) & (subject.times < center + half_width)
            windowed.append(trial_data[:, mask, :].mean(axis=1))
        y = np.stack(windowed, axis=1)
        x, columns = design_matrix(meta)
        beta = np.linalg.pinv(x) @ y.reshape(len(meta), -1)
        beta = beta.reshape(len(columns), len(times), len(channels))
        beta_map = {name: beta[idx] for idx, name in enumerate(columns)}
        maps = {factor: beta_map[factor] for factor in FACTORS}
        maps["Surface"] = np.mean([beta_map[factor] for factor in SURFACE], axis=0)
        maps["Structural"] = np.mean([beta_map[factor] for factor in STRUCTURAL], axis=0)
        maps["Surface_minus_Structural"] = maps["Surface"] - maps["Structural"]
        subject_arrays.append(np.stack([maps[name] for name in CONTRASTS], axis=0))
        used_subjects.append(subject.subj)
        qc["included"] = True
        qc_rows.append(qc)
        logging.info("Computed whole-scalp beta tensor for %s", subject.subj)
    return np.stack(subject_arrays), used_subjects, channels, pd.DataFrame(qc_rows)


def rm_anova_oneway(values: np.ndarray) -> tuple[float, int, int, float]:
    n, k = values.shape
    grand = values.mean()
    subject_means = values.mean(axis=1, keepdims=True)
    stage_means = values.mean(axis=0, keepdims=True)
    ss_stage = n * np.sum((stage_means - grand) ** 2)
    residual = values - subject_means - stage_means + grand
    ss_error = np.sum(residual**2)
    df1 = k - 1
    df2 = (n - 1) * (k - 1)
    f_value = float((ss_stage / df1) / (ss_error / df2))
    return f_value, df1, df2, float(stats.f.sf(f_value, df1, df2))


def rm_permutation_p(values: np.ndarray, n_perm: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    observed = rm_anova_oneway(values)[0]
    null = np.zeros(n_perm)
    for idx in range(n_perm):
        permuted = np.stack([row[rng.permutation(values.shape[1])] for row in values])
        null[idx] = rm_anova_oneway(permuted)[0]
    return float((np.sum(null >= observed) + 1) / (n_perm + 1))


def paired_summary(x: np.ndarray, y: np.ndarray | None = None) -> dict[str, float]:
    values = np.asarray(x, float) if y is None else np.asarray(x, float) - np.asarray(y, float)
    values = values[np.isfinite(values)]
    t, p = stats.ttest_1samp(values, 0.0)
    sd = values.std(ddof=1)
    return {
        "n_subjects": int(len(values)),
        "mean_difference": float(values.mean()),
        "sem_difference": float(sd / math.sqrt(len(values))),
        "t": float(t),
        "p": float(p),
        "cohen_dz": float(values.mean() / sd) if sd > 0 else np.nan,
    }


def stage_interaction(stage_betas: pd.DataFrame, n_perm: int, seed: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    pivot = stage_betas[stage_betas["contrast"] == "Surface_minus_Structural"].pivot(
        index="subj", columns="stage", values="mean_beta_uV"
    )
    pivot = pivot[STAGE_ORDER].dropna()
    values = pivot.to_numpy(float)
    f_value, df1, df2, p_param = rm_anova_oneway(values)
    omnibus = pd.DataFrame(
        [
            {
                "test": "stage_by_edit_class_interaction",
                "n_subjects": len(values),
                "F": f_value,
                "df1": df1,
                "df2": df2,
                "p_parametric": p_param,
                "p_permutation": rm_permutation_p(values, n_perm, seed),
            }
        ]
    )
    rows = []
    for i, stage_a in enumerate(STAGE_ORDER):
        for stage_b in STAGE_ORDER[i + 1 :]:
            summary = paired_summary(pivot[stage_a].to_numpy(), pivot[stage_b].to_numpy())
            rows.append({"stage_a": stage_a, "stage_b": stage_b, **summary})
    pairwise = pd.DataFrame(rows)
    pairwise["p_holm"] = holm(pairwise["p"])
    return omnibus, pairwise


def behavior_interactions(behavior_betas: pd.DataFrame, n_perm: int, seed: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    wide = behavior_betas.pivot_table(
        index="subj", columns=["outcome", "contrast"], values="behavior_beta"
    ).dropna()
    interaction = (
        wide[("Beauty", "Structural")]
        - wide[("Naturalness", "Structural")]
        - wide[("Beauty", "Surface")]
        + wide[("Naturalness", "Surface")]
    )
    rng = np.random.default_rng(seed)
    observed = abs(float(interaction.mean()))
    null = np.array(
        [
            abs(float((interaction.to_numpy() * rng.choice([-1, 1], size=len(interaction))).mean()))
            for _ in range(n_perm)
        ]
    )
    primary = pd.DataFrame(
        [
            {
                "test": "rating_dimension_by_edit_class_interaction",
                **paired_summary(interaction.to_numpy()),
                "p_permutation": float((np.sum(null >= observed) + 1) / (n_perm + 1)),
            }
        ]
    )

    factor_differences = np.column_stack(
        [wide[("Beauty", factor)] - wide[("Naturalness", factor)] for factor in FACTORS]
    )
    f_value, df1, df2, p_param = rm_anova_oneway(factor_differences)
    rows = [
        {
            "test": "rating_dimension_by_specific_factor_interaction",
            "n_subjects": len(wide),
            "F": f_value,
            "df1": df1,
            "df2": df2,
            "p_parametric": p_param,
            "p_permutation": rm_permutation_p(factor_differences, n_perm, seed + 1),
        }
    ]
    for idx, factor in enumerate(FACTORS):
        summary = paired_summary(factor_differences[:, idx])
        rows.append({"test": f"Beauty_minus_Naturalness:{factor}", **summary})
    factor_table = pd.DataFrame(rows)
    mask = factor_table["test"].str.startswith("Beauty_minus")
    factor_table.loc[mask, "p_holm"] = holm(factor_table.loc[mask, "p"])
    return primary, factor_table


def temporal_centroids(beta_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (subj, contrast), sub in beta_df[beta_df["contrast"].isin(["Surface", "Structural"])].groupby(
        ["subj", "contrast"]
    ):
        sub = sub.sort_values("time_ms")
        sub = sub[sub["time_ms"].between(200, 1000)]
        weights = np.abs(sub["beta_uV"].to_numpy(float))
        times = sub["time_ms"].to_numpy(float)
        centroid = float(np.sum(times * weights) / np.sum(weights)) if np.sum(weights) else np.nan
        rows.append({"subj": subj, "contrast": contrast, "temporal_centroid_ms": centroid})
    return pd.DataFrame(rows)


def centroid_test(centroids: pd.DataFrame) -> pd.DataFrame:
    wide = centroids.pivot(index="subj", columns="contrast", values="temporal_centroid_ms").dropna()
    return pd.DataFrame(
        [
            {
                "test": "Structural_minus_Surface_temporal_centroid",
                "surface_mean_ms": float(wide["Surface"].mean()),
                "structural_mean_ms": float(wide["Structural"].mean()),
                **paired_summary(wide["Structural"].to_numpy(), wide["Surface"].to_numpy()),
            }
        ]
    )


def loso_stability(beta_df: pd.DataFrame, clusters: pd.DataFrame) -> pd.DataFrame:
    rows = []
    selected = clusters[clusters["cluster_p"] < 0.05]
    for _, cluster in selected.iterrows():
        contrast = cluster["contrast"]
        start = float(cluster["cluster_start_ms"])
        end = float(cluster["cluster_end_ms"])
        window = (
            beta_df[
                (beta_df["contrast"] == contrast)
                & beta_df["time_ms"].between(start, end)
            ]
            .groupby("subj")["beta_uV"]
            .mean()
        )
        full_sign = np.sign(window.mean())
        pass_count = 0
        signs = 0
        p_values = []
        for held_out in window.index:
            values = window.drop(held_out).to_numpy(float)
            t, p = stats.ttest_1samp(values, 0.0)
            p_values.append(float(p))
            signs += int(np.sign(values.mean()) == full_sign)
            pass_count += int(p < 0.05 and np.sign(values.mean()) == full_sign)
        rows.append(
            {
                "contrast": contrast,
                "window": f"{start:.0f}-{end:.0f} ms",
                "n_leave_one_out_runs": int(len(window)),
                "same_direction_fraction": float(signs / len(window)),
                "p_lt_0.05_same_direction_fraction": float(pass_count / len(window)),
                "median_loso_p": float(np.median(p_values)),
                "max_loso_p": float(np.max(p_values)),
            }
        )
    return pd.DataFrame(rows)


def cluster_1d_matrix(data: np.ndarray, times: np.ndarray, n_perm: int, seed: int) -> list[dict[str, float]]:
    rng = np.random.default_rng(seed)
    t_obs = t_values(data)
    threshold = float(stats.t.ppf(0.975, data.shape[0] - 1))

    def clusters(tvals: np.ndarray) -> list[np.ndarray]:
        out = []
        for mask in (tvals >= threshold, tvals <= -threshold):
            padded = np.r_[False, mask, False]
            starts = np.where(np.diff(padded.astype(int)) == 1)[0]
            ends = np.where(np.diff(padded.astype(int)) == -1)[0]
            out.extend([np.arange(start, end) for start, end in zip(starts, ends)])
        return out

    observed = clusters(t_obs)
    null = np.zeros(n_perm)
    for idx in range(n_perm):
        signs = rng.choice([-1.0, 1.0], size=(data.shape[0], 1))
        perm_t = t_values(data * signs)
        null[idx] = max((np.sum(np.abs(perm_t[c])) for c in clusters(perm_t)), default=0.0)
    rows = []
    for cluster in observed:
        mass = float(np.sum(np.abs(t_obs[cluster])))
        rows.append(
            {
                "cluster_start_ms": float(times[cluster[0]]),
                "cluster_end_ms": float(times[cluster[-1]]),
                "cluster_mass_abs_t": mass,
                "cluster_p": float((np.sum(null >= mass) + 1) / (n_perm + 1)),
                "peak_time_ms": float(times[cluster[np.argmax(np.abs(t_obs[cluster]))]]),
            }
        )
    return rows


def robustness_checks(
    tensor: np.ndarray,
    times: np.ndarray,
    channels: list[str],
    n_perm: int,
    seed: int,
) -> pd.DataFrame:
    roi_defs = {
        "Combined_posterior_CP": [ch for ch in N170_LPP_CHANNELS if ch in channels],
        "Posterior": [ch for ch in POSTERIOR_CHANNELS if ch in channels],
        "Occipitotemporal": [
            ch for ch in ["P7", "P5", "P6", "P8", "PO7", "PO5", "PO6", "PO8", "O1", "O2"] if ch in channels
        ],
        "Centroparietal": [
            ch for ch in ["CP5", "CP3", "CP1", "CPZ", "CP2", "CP4", "CP6", "P3", "P1", "PZ", "P2", "P4"] if ch in channels
        ],
    }
    rows = []
    contrast_lookup = {contrast: idx for idx, contrast in enumerate(CONTRASTS)}
    step = float(np.median(np.diff(times)))
    for roi_idx, (roi_name, roi_channels) in enumerate(roi_defs.items()):
        indices = [channels.index(ch) for ch in roi_channels]
        for smooth_ms in (20, 40, 60):
            width = max(1, int(round(smooth_ms / step)))
            for contrast in PLOT_CONTRASTS:
                data = tensor[:, contrast_lookup[contrast], :, :][:, :, indices].mean(axis=2)
                data = moving_average(data, width)
                cluster_rows = cluster_1d_matrix(
                    data,
                    times,
                    n_perm=n_perm,
                    seed=seed + roi_idx * 100 + smooth_ms + contrast_lookup[contrast],
                )
                best = min(cluster_rows, key=lambda item: item["cluster_p"], default=None)
                rows.append(
                    {
                        "roi": roi_name,
                        "n_channels": len(indices),
                        "smoothing_ms": smooth_ms,
                        "contrast": contrast,
                        "best_cluster_start_ms": best["cluster_start_ms"] if best else np.nan,
                        "best_cluster_end_ms": best["cluster_end_ms"] if best else np.nan,
                        "best_cluster_p": best["cluster_p"] if best else np.nan,
                        "has_corrected_cluster": bool(best and best["cluster_p"] < 0.05),
                    }
                )
    return pd.DataFrame(rows)


def plot_spatiotemporal_maps(
    t_maps: dict[str, np.ndarray],
    sig_maps: dict[str, np.ndarray],
    times: np.ndarray,
    channels: list[str],
    outdir: Path,
) -> None:
    selected = ["Surface", "Structural", "Surface_minus_Structural"]
    vmax = max(3.0, max(float(np.nanmax(np.abs(t_maps[name]))) for name in selected))
    fig, axes = plt.subplots(len(selected), 1, figsize=(13, 9), sharex=True)
    for ax, name in zip(axes, selected):
        image = ax.imshow(
            t_maps[name].T,
            aspect="auto",
            origin="lower",
            extent=[times.min(), times.max(), -0.5, len(channels) - 0.5],
            cmap="RdBu_r",
            vmin=-vmax,
            vmax=vmax,
        )
        sig = sig_maps[name].T > 0
        if sig.any():
            ax.contour(times, np.arange(len(channels)), sig, levels=[0.5], colors="black", linewidths=0.8)
        tick_idx = np.linspace(0, len(channels) - 1, min(12, len(channels)), dtype=int)
        ax.set_yticks(tick_idx, [channels[idx] for idx in tick_idx], fontsize=8)
        ax.set_ylabel(name.replace("_", "\n"))
        for boundary in (140, 200, 300, 450, 800):
            ax.axvline(boundary, color="#111827", lw=0.5, alpha=0.35)
    axes[-1].set_xlabel("Time from face onset (ms)")
    fig.colorbar(image, ax=axes, shrink=0.75, label="one-sample t")
    fig.suptitle("Whole-scalp spatiotemporal factor effects\nblack contours: cluster-corrected p < .05", y=0.995)
    save_figure(fig, outdir / "Fig2_whole_scalp_spatiotemporal_clusters")


def plot_stage_interaction(stage_betas: pd.DataFrame, omnibus: pd.DataFrame, outdir: Path) -> None:
    focus = stage_betas[stage_betas["contrast"].isin(["Surface", "Structural"])].copy()
    summary = (
        focus.groupby(["stage", "contrast"])["mean_beta_uV"]
        .agg(["mean", "sem"])
        .reset_index()
    )
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    x = np.arange(len(STAGE_ORDER))
    for contrast, offset in (("Surface", -0.08), ("Structural", 0.08)):
        sub = summary[summary["contrast"] == contrast].set_index("stage").reindex(STAGE_ORDER)
        ax.errorbar(
            x + offset,
            sub["mean"],
            yerr=sub["sem"],
            marker="o",
            ms=7,
            lw=2,
            capsize=4,
            color=COLORS[contrast],
            label=contrast,
        )
    row = omnibus.iloc[0]
    ax.text(
        0.02,
        0.97,
        f"Stage x edit class: F({int(row.df1)}, {int(row.df2)})={row.F:.2f}, "
        f"permutation p={row.p_permutation:.3g}",
        transform=ax.transAxes,
        ha="left",
        va="top",
    )
    ax.axhline(0, color="black", ls="--", lw=1)
    ax.set_xticks(x, [stage.replace("_", "\n") for stage in STAGE_ORDER])
    ax.set_ylabel("Factor-effect beta (uV)")
    ax.set_title("Neural processing stage differentiates surface and structural edits")
    ax.legend(frameon=False)
    fig.tight_layout()
    save_figure(fig, outdir / "Fig3_stage_by_edit_class_interaction")


def plot_behavior_interaction(
    behavior_betas: pd.DataFrame,
    primary: pd.DataFrame,
    outdir: Path,
) -> None:
    group = (
        behavior_betas[behavior_betas["contrast"].isin(["Surface", "Structural"])]
        .groupby(["outcome", "contrast"])["behavior_beta"]
        .agg(["mean", "sem"])
        .reset_index()
    )
    factor_group = (
        behavior_betas[behavior_betas["contrast"].isin(FACTORS)]
        .groupby(["outcome", "contrast"])["behavior_beta"]
        .mean()
        .unstack("outcome")
        .reindex(FACTORS)
    )
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    x = np.arange(2)
    for contrast, offset in (("Surface", -0.06), ("Structural", 0.06)):
        sub = group[group["contrast"] == contrast].set_index("outcome").reindex(["Naturalness", "Beauty"])
        axes[0].errorbar(
            x + offset,
            sub["mean"],
            yerr=sub["sem"],
            marker="o",
            lw=2,
            ms=7,
            capsize=4,
            label=contrast,
            color=COLORS[contrast],
        )
    p_perm = primary.iloc[0]["p_permutation"]
    axes[0].text(0.02, 0.97, f"Dimension x edit class permutation p={p_perm:.3g}", transform=axes[0].transAxes, va="top")
    axes[0].axhline(0, color="black", ls="--", lw=1)
    axes[0].set_xticks(x, ["Naturalness", "Beauty"])
    axes[0].set_ylabel("Rating beta")
    axes[0].set_title("Beauty-naturalness dissociation")
    axes[0].legend(frameon=False)

    matrix = factor_group[["Naturalness", "Beauty"]].to_numpy(float)
    vmax = max(0.2, float(np.nanmax(np.abs(matrix))))
    image = axes[1].imshow(matrix, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    axes[1].set_xticks([0, 1], ["Naturalness", "Beauty"])
    axes[1].set_yticks(np.arange(len(FACTORS)), FACTORS)
    for row in range(matrix.shape[0]):
        for col in range(matrix.shape[1]):
            axes[1].text(col, row, f"{matrix[row, col]:.2f}", ha="center", va="center", fontsize=10)
    axes[1].set_title("Specific editing-factor effects")
    fig.colorbar(image, ax=axes[1], fraction=0.046, label="rating beta")
    fig.tight_layout()
    save_figure(fig, outdir / "Fig4_behavioral_dissociation")


def plot_robustness(robustness: pd.DataFrame, loso: pd.DataFrame, outdir: Path) -> None:
    contrast_order = PLOT_CONTRASTS
    configs = (
        robustness.assign(config=lambda d: d["roi"] + "\n" + d["smoothing_ms"].astype(str) + " ms")
        .pivot(index="config", columns="contrast", values="best_cluster_p")
    )
    configs = configs.reindex(columns=contrast_order)
    score = -np.log10(configs.to_numpy(float))
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={"width_ratios": [1.7, 1]})
    image = axes[0].imshow(score, aspect="auto", cmap="viridis", vmin=0, vmax=max(2, np.nanmax(score)))
    axes[0].set_xticks(np.arange(len(contrast_order)), [name.replace("_", "\n") for name in contrast_order])
    axes[0].set_yticks(np.arange(len(configs)), configs.index, fontsize=8)
    axes[0].set_title("ROI and smoothing robustness")
    for row in range(score.shape[0]):
        for col in range(score.shape[1]):
            p_value = configs.iloc[row, col]
            axes[0].text(col, row, "n/a" if pd.isna(p_value) else f"{p_value:.3f}", ha="center", va="center", fontsize=7, color="white" if score[row, col] > 1 else "black")
    fig.colorbar(image, ax=axes[0], fraction=0.035, label="-log10(cluster p)")

    if not loso.empty:
        y = np.arange(len(loso))
        axes[1].barh(y, loso["p_lt_0.05_same_direction_fraction"], color=[COLORS.get(c, "#64748b") for c in loso["contrast"]])
        axes[1].set_yticks(y, loso["contrast"].str.replace("_", " "))
        axes[1].set_xlim(0, 1)
        axes[1].set_xlabel("LOSO runs retaining p < .05 and direction")
        axes[1].set_title("Leave-one-subject-out stability")
    fig.tight_layout()
    save_figure(fig, outdir / "Fig5_robustness_and_loso")


def plot_integrated_timeline(
    group_stats: pd.DataFrame,
    clusters: pd.DataFrame,
    outdir: Path,
) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(13, 7), sharex=True, gridspec_kw={"height_ratios": [3, 1.4]})
    for contrast in ["Mouth", "Skin", "Structural", "Surface"]:
        sub = group_stats[group_stats["contrast"] == contrast]
        axes[0].plot(sub["time_ms"], sub["mean_beta_uV"], lw=2, color=COLORS[contrast], label=contrast)
        axes[0].fill_between(sub["time_ms"], sub["ci95_low"], sub["ci95_high"], color=COLORS[contrast], alpha=0.12)
    axes[0].axhline(0, color="black", ls="--", lw=1)
    axes[0].set_ylabel("Factor-effect beta (uV)")
    axes[0].legend(frameon=False, ncol=4)
    axes[0].set_title("Integrated temporal account of facial-edit processing")
    significant = clusters[clusters["cluster_p"] < 0.10].copy()
    y_order = ["Mouth", "Skin", "Surface", "Surface_minus_Structural", "Structural"]
    for y_idx, contrast in enumerate(y_order):
        subset = significant[significant["contrast"] == contrast]
        for _, row in subset.iterrows():
            axes[1].hlines(y_idx, row["cluster_start_ms"], row["cluster_end_ms"], color=COLORS[contrast], lw=9, alpha=0.95 if row["cluster_p"] < 0.05 else 0.45)
            axes[1].text(row["cluster_end_ms"] + 8, y_idx, f"p={row['cluster_p']:.3f}", va="center", fontsize=8)
    axes[1].set_yticks(np.arange(len(y_order)), [name.replace("_", " ") for name in y_order])
    axes[1].set_xlabel("Time from face onset (ms)")
    axes[1].set_xlim(0, 1050)
    axes[1].grid(axis="x", alpha=0.2)
    fig.tight_layout()
    save_figure(fig, outdir / "Fig1_integrated_temporal_results")


def dataframe_records(df: pd.DataFrame) -> list[dict[str, Any]]:
    return json.loads(df.to_json(orient="records"))


def main() -> None:
    project_root = Path(os.environ.get("IPM_DATA_ROOT", Path.cwd())).resolve()
    out = project_root / "paper_extension_final"
    dirs = setup_dirs(out)
    setup_logging(dirs["logs"] / "run_paper_grade_extension.log")
    seed = 20260731
    n_perm_spatiotemporal = 2000
    n_perm_other = 5000
    logging.info("Starting paper-grade extension in %s", out)

    eeg_dir = project_root / "erp_dynamics_final" / "tables"
    behavior_dir = project_root / "erp_behavior_integrated_final" / "tables"
    beta_df = pd.read_csv(eeg_dir / "erp_subject_time_resolved_betas.csv")
    stage_betas = pd.read_csv(eeg_dir / "erp_subject_stage_betas.csv")
    group_stats = pd.read_csv(eeg_dir / "erp_time_resolved_group_stats.csv")
    temporal_clusters = pd.read_csv(eeg_dir / "erp_time_resolved_clusters.csv")
    behavior_betas = pd.read_csv(behavior_dir / "behavior_subject_factor_betas.csv")

    centers = np.arange(0, 1001, 20, dtype=float)
    tensor_path = dirs["intermediate"] / "whole_scalp_subject_betas.npz"
    if tensor_path.exists():
        loaded = np.load(tensor_path, allow_pickle=True)
        tensor = loaded["tensor"]
        eeg_subjects = loaded["subjects"].astype(str).tolist()
        channels = loaded["channels"].astype(str).tolist()
        qc = pd.read_csv(dirs["tables"] / "whole_scalp_data_qc.csv")
        logging.info("Loaded cached whole-scalp tensor %s", tensor.shape)
    else:
        tensor, eeg_subjects, channels, qc = load_channel_betas(project_root, centers)
        np.savez_compressed(
            tensor_path,
            tensor=tensor.astype(np.float32),
            subjects=np.asarray(eeg_subjects),
            channels=np.asarray(channels),
            times=centers,
            contrasts=np.asarray(CONTRASTS),
        )
        save_csv(qc, dirs["tables"] / "whole_scalp_data_qc.csv")

    adjacency = build_adjacency(channels)
    cluster_parts = []
    t_maps: dict[str, np.ndarray] = {}
    sig_maps: dict[str, np.ndarray] = {}
    for contrast_idx, contrast in enumerate(CONTRASTS):
        logging.info("Spatiotemporal cluster test: %s", contrast)
        cluster_df, t_map, sig_map = spatiotemporal_cluster_test(
            tensor[:, contrast_idx],
            centers,
            channels,
            adjacency,
            n_perm=n_perm_spatiotemporal,
            seed=seed + contrast_idx,
        )
        if not cluster_df.empty:
            cluster_df.insert(0, "contrast", contrast)
            cluster_parts.append(cluster_df)
        t_maps[contrast] = t_map
        sig_maps[contrast] = sig_map
    spatial_clusters = pd.concat(cluster_parts, ignore_index=True) if cluster_parts else pd.DataFrame()
    if not spatial_clusters.empty:
        min_p = spatial_clusters.groupby("contrast")["cluster_p"].min()
        adjusted = dict(zip(min_p.index, holm(min_p.to_numpy())))
        spatial_clusters["min_cluster_p_holm_across_contrasts"] = spatial_clusters["contrast"].map(adjusted)
        spatial_clusters = spatial_clusters.sort_values(["cluster_p", "contrast"]).reset_index(drop=True)
    save_csv(spatial_clusters, dirs["tables"] / "whole_scalp_spatiotemporal_clusters.csv")

    omnibus, stage_pairwise = stage_interaction(stage_betas, n_perm_other, seed)
    save_csv(omnibus, dirs["tables"] / "stage_by_edit_class_interaction.csv")
    save_csv(stage_pairwise, dirs["tables"] / "stage_interaction_pairwise.csv")

    behavior_primary, behavior_factor = behavior_interactions(behavior_betas, n_perm_other, seed)
    save_csv(behavior_primary, dirs["tables"] / "behavior_dimension_by_edit_class_interaction.csv")
    save_csv(behavior_factor, dirs["tables"] / "behavior_dimension_by_factor_interactions.csv")
    matched_behavior = behavior_betas[behavior_betas["subj"].isin(eeg_subjects)].copy()
    matched_primary, matched_factor = behavior_interactions(matched_behavior, n_perm_other, seed + 100)
    save_csv(
        matched_primary,
        dirs["tables"] / "behavior_dimension_by_edit_class_interaction_eeg_matched.csv",
    )
    save_csv(
        matched_factor,
        dirs["tables"] / "behavior_dimension_by_factor_interactions_eeg_matched.csv",
    )

    centroids = temporal_centroids(beta_df)
    centroid_stats = centroid_test(centroids)
    save_csv(centroids, dirs["tables"] / "subject_temporal_centroids.csv")
    save_csv(centroid_stats, dirs["tables"] / "temporal_centroid_comparison.csv")

    loso = loso_stability(beta_df, temporal_clusters)
    save_csv(loso, dirs["tables"] / "leave_one_subject_out_stability.csv")

    robustness = robustness_checks(
        tensor,
        centers,
        channels,
        n_perm=500,
        seed=seed,
    )
    save_csv(robustness, dirs["tables"] / "roi_smoothing_robustness.csv")

    plot_integrated_timeline(group_stats, temporal_clusters, dirs["figures"])
    plot_spatiotemporal_maps(t_maps, sig_maps, centers, channels, dirs["figures"])
    plot_stage_interaction(stage_betas, omnibus, dirs["figures"])
    plot_behavior_interaction(behavior_betas, behavior_primary, dirs["figures"])
    plot_robustness(robustness, loso, dirs["figures"])

    summary = {
        "analysis": "paper-grade extension for edited-face EEG",
        "seed": seed,
        "n_eeg_subjects_whole_scalp": len(eeg_subjects),
        "eeg_subjects_whole_scalp": eeg_subjects,
        "n_behavior_subjects": int(behavior_betas["subj"].nunique()),
        "n_channels_whole_scalp": len(channels),
        "whole_scalp_clusters_p_lt_0_05": dataframe_records(
            spatial_clusters[spatial_clusters["cluster_p"] < 0.05]
        )
        if not spatial_clusters.empty
        else [],
        "stage_interaction": dataframe_records(omnibus),
        "behavior_interaction": dataframe_records(behavior_primary),
        "behavior_interaction_eeg_matched": dataframe_records(matched_primary),
        "temporal_centroid": dataframe_records(centroid_stats),
        "loso_stability": dataframe_records(loso),
        "robustness_significant_config_counts": dataframe_records(
            robustness.groupby("contrast")["has_corrected_cluster"].agg(["sum", "count"]).reset_index()
        ),
        "figures": [str(path.resolve()) for path in sorted(dirs["figures"].glob("*.png"))],
        "limitations": [
            "The whole-scalp extension is confirmatory for spatial distribution but does not provide source localization.",
            "Cross-subject EEG-behavior correlations remain exploratory and were not strengthened by this extension.",
            "Temporal cluster onset and duration describe corrected evidence, not the exact onset or duration of a latent cognitive process.",
        ],
    }
    summary_path = dirs["summaries"] / "paper_extension_results.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, default=json_default), encoding="utf-8")
    logging.info("Finished paper-grade extension. Summary: %s", summary_path)


if __name__ == "__main__":
    main()
