from pathlib import Path
import re
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.size": 11,
    "axes.labelsize": 12,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
})

from scipy.stats import fisher_exact, mannwhitneyu
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.proportion import proportion_confint

warnings.filterwarnings("ignore")

base = Path(__file__).resolve().parents[1]
meta_dir = base / "data/raw/metadata/metadata"
prof_dir = base / "data/raw/metaphlan4_profiles/metaphlan4_profiles"
res_dir = base / "results"
fig_dir = base / "figures"
res_dir.mkdir(exist_ok=True, parents=True)
fig_dir.mkdir(exist_ok=True, parents=True)

studies = [
    "Public_study__FengQ_2015",
    "Public_study__GuptaA_2019",
    "Public_study__VogtmannE_2016",
    "Public_study__WirbelJ_2018",
    "Public_study__YachidaS_2019",
    "Public_study__YangJ_2020",
    "Public_study__YuJ_2015",
    "Public_study__ZellerG_2014",
]

GM = "Gemella morbillorum"
GROUP_ORDER = ["Healthy", "Stage 0", "Stage I", "Stage II", "Stage III", "Stage IV"]


def clean_species_name(x):
    x = str(x)
    if "|" in x:
        x = x.split("|")[-1]
    x = x.replace("s__", "")
    x = x.replace("_", " ")
    return x


def read_metaphlan_species_table(path):
    df = pd.read_csv(path, sep="\t", comment="#")
    first_col = df.columns[0]
    df = df.rename(columns={first_col: "clade_name"})
    df = df[df["clade_name"].astype(str).str.contains("s__")]
    df = df[~df["clade_name"].astype(str).str.contains("t__")]
    df["species"] = df["clade_name"].map(clean_species_name)
    df = df.drop(columns=["clade_name"])
    df = df.set_index("species")
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0)
    df = df.T
    df.index.name = "sample_id"
    return df


def normalize_stage(x):
    if pd.isna(x):
        return np.nan

    s = str(x).strip().lower()
    s = s.replace("stage", "")
    s = s.replace("ajcc", "")
    s = s.replace(" ", "")
    s = s.replace("-", "")
    s = s.replace("_", "")

    if s in ["0", "zero"]:
        return "Stage 0"

    if re.search(r"\biv\b|^iv|4", s):
        return "Stage IV"
    if re.search(r"\biii\b|^iii|3", s):
        return "Stage III"
    if re.search(r"\bii\b|^ii|2", s):
        return "Stage II"
    if re.search(r"\bi\b|^i|1", s):
        return "Stage I"

    return np.nan


def call_group(row):
    text = " ".join([str(v).lower() for v in row.values])
    healthy_terms = ["healthy", "control", "normal"]
    crc_terms = ["colorectal cancer", "crc", "colon cancer", "rectal cancer", "adenocarcinoma"]

    if any(t in text for t in healthy_terms) and not any(t in text for t in crc_terms):
        return "Healthy"

    if any(t in text for t in crc_terms):
        return row.get("stage_group", np.nan)

    return np.nan


def get_sample_col(meta):
    candidates = [
        "Name",
        "Sample Source ID",
        "sample_id",
        "sample",
        "Sample",
        "sample_alias",
    ]
    for c in candidates:
        if c in meta.columns:
            return c
    return meta.columns[0]


all_rows = []

for study in studies:
    meta_path = meta_dir / f"{study}.tsv"
    prof_path = prof_dir / f"{study}.tsv"

    meta = pd.read_csv(meta_path, sep="\t")
    prof = read_metaphlan_species_table(prof_path)

    sample_col = get_sample_col(meta)
    meta = meta.copy()
    meta["sample_id"] = meta[sample_col].astype(str)
    meta["study"] = study.replace("Public_study__", "")

    stage_col = None
    for c in ["Tumor Staging AJCC", "tumor_stage", "stage", "Stage"]:
        if c in meta.columns:
            stage_col = c
            break

    if stage_col is None:
        meta["stage_group"] = np.nan
    else:
        meta["stage_group"] = meta[stage_col].map(normalize_stage)

    meta["analysis_group"] = meta.apply(call_group, axis=1)

    merged = meta.merge(
        prof,
        left_on="sample_id",
        right_index=True,
        how="inner",
    )

    all_rows.append(merged)

dat = pd.concat(all_rows, axis=0, ignore_index=True)
dat = dat[dat["analysis_group"].isin(GROUP_ORDER)].copy()

if GM not in dat.columns:
    raise ValueError(f"{GM} was not found in the species table.")

dat["gm_abundance"] = pd.to_numeric(dat[GM], errors="coerce").fillna(0)
dat["gm_detected"] = dat["gm_abundance"] > 0

rng = np.random.default_rng(42)


def bootstrap_mean_ci(x, n_boot=5000, alpha=0.05):
    x = np.asarray(x, dtype=float)
    if len(x) == 0:
        return np.nan, np.nan
    means = []
    for _ in range(n_boot):
        sample = rng.choice(x, size=len(x), replace=True)
        means.append(sample.mean())
    return (
        np.quantile(means, alpha / 2),
        np.quantile(means, 1 - alpha / 2),
    )


summary_rows = []

for group in GROUP_ORDER:
    sub = dat[dat["analysis_group"] == group]
    n = len(sub)
    n_det = int(sub["gm_detected"].sum())
    detection_rate = n_det / n if n else np.nan

    ci_low, ci_high = proportion_confint(
        count=n_det,
        nobs=n,
        alpha=0.05,
        method="wilson",
    )

    mean_abundance = sub["gm_abundance"].mean()
    mean_low, mean_high = bootstrap_mean_ci(sub["gm_abundance"].values)

    summary_rows.append({
        "group": group,
        "n": n,
        "n_detected": n_det,
        "detection_rate": detection_rate,
        "detection_rate_percent": detection_rate * 100,
        "detection_ci_low": ci_low,
        "detection_ci_high": ci_high,
        "detection_ci_low_percent": ci_low * 100,
        "detection_ci_high_percent": ci_high * 100,
        "mean_abundance": mean_abundance,
        "mean_abundance_ci_low": mean_low,
        "mean_abundance_ci_high": mean_high,
        "median_abundance": sub["gm_abundance"].median(),
        "positive_only_mean_abundance": sub.loc[sub["gm_detected"], "gm_abundance"].mean(),
    })

summary = pd.DataFrame(summary_rows)
summary.to_csv(res_dir / "02_gm_stage_summary_with_errorbars.csv", index=False)

pair_rows = []

for i, g1 in enumerate(GROUP_ORDER):
    for g2 in GROUP_ORDER[i + 1:]:
        a = dat[dat["analysis_group"] == g1]
        b = dat[dat["analysis_group"] == g2]

        table = [
            [int(a["gm_detected"].sum()), int((~a["gm_detected"]).sum())],
            [int(b["gm_detected"].sum()), int((~b["gm_detected"]).sum())],
        ]

        try:
            odds_ratio, p_det = fisher_exact(table)
        except ValueError:
            odds_ratio, p_det = np.nan, 1.0

        try:
            stat_abun, p_abun = mannwhitneyu(
                a["gm_abundance"],
                b["gm_abundance"],
                alternative="two-sided",
            )
        except ValueError:
            stat_abun, p_abun = np.nan, 1.0

        pair_rows.append({
            "group1": g1,
            "group2": g2,
            "test_detection": "Fisher exact test",
            "odds_ratio_group1_vs_group2": odds_ratio,
            "p_detection": p_det,
            "test_abundance": "Mann-Whitney U test",
            "mannwhitney_u": stat_abun,
            "p_abundance": p_abun,
            "detection_rate_group1": a["gm_detected"].mean(),
            "detection_rate_group2": b["gm_detected"].mean(),
            "mean_abundance_group1": a["gm_abundance"].mean(),
            "mean_abundance_group2": b["gm_abundance"].mean(),
        })

pairs = pd.DataFrame(pair_rows)
pairs["fdr_detection"] = multipletests(pairs["p_detection"], method="fdr_bh")[1]
pairs["fdr_abundance"] = multipletests(pairs["p_abundance"], method="fdr_bh")[1]

pairs.to_csv(res_dir / "02_gm_stage_pairwise_tests.csv", index=False)
pairs[[
    "group1", "group2", "odds_ratio_group1_vs_group2",
    "p_detection", "fdr_detection",
    "detection_rate_group1", "detection_rate_group2",
]].to_csv(res_dir / "02_gm_stage_pairwise_detection_fisher.csv", index=False)

pairs[[
    "group1", "group2", "mannwhitney_u",
    "p_abundance", "fdr_abundance",
    "mean_abundance_group1", "mean_abundance_group2",
]].to_csv(res_dir / "02_gm_stage_pairwise_abundance_mannwhitney.csv", index=False)


def plot_detection(summary):
    d = summary.copy()
    x = np.arange(len(d))
    y = d["detection_rate_percent"].values
    yerr = np.vstack([
        y - d["detection_ci_low_percent"].values,
        d["detection_ci_high_percent"].values - y,
    ])

    plt.figure(figsize=(7, 4.5))
    plt.bar(x, y, yerr=yerr, capsize=4)
    plt.xticks(x, d["group"], rotation=30, ha="right")
    plt.ylabel("Detection rate of G. morbillorum (%)")
    plt.xlabel("")
    plt.tight_layout()
    plt.savefig(fig_dir / "gm_detection_rate_by_healthy_stage.png", dpi=300)
    plt.savefig(fig_dir / "gm_detection_rate_by_healthy_stage.pdf")
    plt.close()


def plot_abundance(summary):
    d = summary.copy()
    x = np.arange(len(d))
    y = d["mean_abundance"].values
    yerr = np.vstack([
        y - d["mean_abundance_ci_low"].values,
        d["mean_abundance_ci_high"].values - y,
    ])

    plt.figure(figsize=(7, 4.5))
    plt.bar(x, y, yerr=yerr, capsize=4)
    plt.xticks(x, d["group"], rotation=30, ha="right")
    plt.ylabel("Mean relative abundance of G. morbillorum")
    plt.xlabel("")
    plt.tight_layout()
    plt.savefig(fig_dir / "gm_mean_abundance_by_healthy_stage.png", dpi=300)
    plt.savefig(fig_dir / "gm_mean_abundance_by_healthy_stage.pdf")
    plt.close()


plot_detection(summary)
plot_abundance(summary)

print("[INFO] G. morbillorum stage summary")
print(summary.to_string(index=False))
print("[INFO] Done: G. morbillorum stage detection and abundance")

# BEGIN_SIGNIFICANCE_ANNOTATED_PLOTS
# Add FDR-based significance labels for Healthy vs each CRC stage.
# Detection rate: Fisher's exact test with Benjamini-Hochberg FDR correction.
# Mean abundance: Mann-Whitney U test with Benjamini-Hochberg FDR correction.

def _sig_label(fdr):
    if pd.isna(fdr):
        return "NA"
    if fdr < 0.0001:
        return "****"
    if fdr < 0.001:
        return "***"
    if fdr < 0.01:
        return "**"
    if fdr < 0.05:
        return "*"
    return "ns"


def _get_vs_healthy_fdr(pairwise_df, stage, fdr_col):
    rows = pairwise_df[
        ((pairwise_df["group1"] == "Healthy") & (pairwise_df["group2"] == stage))
        | ((pairwise_df["group2"] == "Healthy") & (pairwise_df["group1"] == stage))
    ]
    if rows.empty:
        return np.nan
    return float(rows.iloc[0][fdr_col])


def _plot_detection_with_significance():
    base_dir = Path(__file__).resolve().parents[1]
    results_dir = base_dir / "results"
    figures_dir = base_dir / "figures"

    summary = pd.read_csv(results_dir / "02_gm_stage_summary_with_errorbars.csv")
    pairwise = pd.read_csv(results_dir / "02_gm_stage_pairwise_tests.csv")

    group_order = ["Healthy", "Stage 0", "Stage I", "Stage II", "Stage III", "Stage IV"]
    summary = summary.set_index("group").loc[group_order].reset_index()

    x = np.arange(len(summary))
    y = summary["detection_rate_percent"].values
    yerr_lower = y - summary["detection_ci_low_percent"].values
    yerr_upper = summary["detection_ci_high_percent"].values - y

    fig, ax = plt.subplots(figsize=(7, 4.8))
    ax.bar(
        x, y,
        yerr=[np.zeros_like(yerr_upper), yerr_upper],
        capsize=4,
        color="#8FA1BA",
        edgecolor="#3F5D85",
        ecolor="#3F5D85",
        linewidth=1.2,
        error_kw={"elinewidth": 1.2, "capthick": 1.2}
    )

    ymax = max(y + yerr_upper)
    offset = ymax * 0.06 if ymax > 0 else 1

    for i, group in enumerate(summary["group"]):
        if group == "Healthy":
            continue
        fdr = _get_vs_healthy_fdr(pairwise, group, "fdr_detection")
        label = _sig_label(fdr)
        ax.text(
            i,
            y[i] + yerr_upper[i] + offset * 0.5,
            label,
            ha="center",
            va="bottom",
            fontsize=12,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(summary["group"], rotation=30, ha="right")
    ax.set_ylabel("Detection rate of G. morbillorum (%)")
    ax.set_xlabel("")
    ax.set_ylim(0, max(y + yerr_upper) + offset * 4)


    fig.tight_layout()
    fig.savefig(figures_dir / "gm_detection_rate_by_healthy_stage.png", dpi=300)
    fig.savefig(figures_dir / "gm_detection_rate_by_healthy_stage.pdf")
    plt.close(fig)


def _plot_abundance_with_significance():
    base_dir = Path(__file__).resolve().parents[1]
    results_dir = base_dir / "results"
    figures_dir = base_dir / "figures"

    summary = pd.read_csv(results_dir / "02_gm_stage_summary_with_errorbars.csv")
    pairwise = pd.read_csv(results_dir / "02_gm_stage_pairwise_tests.csv")

    group_order = ["Healthy", "Stage 0", "Stage I", "Stage II", "Stage III", "Stage IV"]
    summary = summary.set_index("group").loc[group_order].reset_index()

    x = np.arange(len(summary))
    y = summary["mean_abundance"].values
    yerr_lower = y - summary["mean_abundance_ci_low"].values
    yerr_upper = summary["mean_abundance_ci_high"].values - y

    fig, ax = plt.subplots(figsize=(7, 4.8))
    ax.bar(
        x, y,
        yerr=[np.zeros_like(yerr_upper), yerr_upper],
        capsize=4,
        color="#8FA1BA",
        edgecolor="#3F5D85",
        ecolor="#3F5D85",
        linewidth=1.2,
        error_kw={"elinewidth": 1.2, "capthick": 1.2}
    )

    ymax = max(y + yerr_upper)
    offset = ymax * 0.08 if ymax > 0 else 0.001

    for i, group in enumerate(summary["group"]):
        if group == "Healthy":
            continue
        fdr = _get_vs_healthy_fdr(pairwise, group, "fdr_abundance")
        label = _sig_label(fdr)
        ax.text(
            i,
            y[i] + yerr_upper[i] + offset * 0.5,
            label,
            ha="center",
            va="bottom",
            fontsize=12,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(summary["group"], rotation=30, ha="right")
    ax.set_ylabel("Mean relative abundance of G. morbillorum")
    ax.set_xlabel("")
    ax.set_ylim(0, max(y + yerr_upper) + offset * 4)


    fig.tight_layout()
    fig.savefig(figures_dir / "gm_mean_abundance_by_healthy_stage.png", dpi=300)
    fig.savefig(figures_dir / "gm_mean_abundance_by_healthy_stage.pdf")
    plt.close(fig)


_plot_detection_with_significance()
_plot_abundance_with_significance()
print("[INFO] Updated G. morbillorum plots with FDR-based significance labels.")
# END_SIGNIFICANCE_ANNOTATED_PLOTS
