from pathlib import Path
import re
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

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


def is_crc_row(row):
    text = " ".join([str(v).lower() for v in row.values])
    crc_terms = [
        "colorectal cancer",
        "crc",
        "colon cancer",
        "rectal cancer",
        "adenocarcinoma",
    ]
    healthy_terms = [
        "healthy",
        "control",
        "normal",
    ]
    if any(t in text for t in crc_terms):
        return True
    if any(t in text for t in healthy_terms):
        return False
    return False


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

    meta["is_crc"] = meta.apply(is_crc_row, axis=1)

    merged = meta.merge(
        prof,
        left_on="sample_id",
        right_index=True,
        how="inner",
    )

    all_rows.append(merged)

dat = pd.concat(all_rows, axis=0, ignore_index=True)

species_cols = [
    c for c in dat.columns
    if c not in set(pd.concat([
        pd.read_csv(meta_dir / f"{s}.tsv", sep="\t").head(0)
        for s in studies
    ], axis=0).columns)
    and c not in ["sample_id", "study", "stage_group", "is_crc"]
]

species_cols = [c for c in species_cols if pd.api.types.is_numeric_dtype(dat[c])]

crc = dat[(dat["is_crc"]) & (dat["stage_group"].notna())].copy()

crc["stage_binary"] = np.nan
crc.loc[crc["stage_group"].isin(["Stage 0", "Stage I"]), "stage_binary"] = 0
crc.loc[crc["stage_group"].isin(["Stage II", "Stage III", "Stage IV"]), "stage_binary"] = 1
crc = crc[crc["stage_binary"].notna()].copy()
crc["stage_binary"] = crc["stage_binary"].astype(int)

X_raw = crc[species_cols].copy()
detected = (X_raw > 0).sum(axis=0)
var = X_raw.var(axis=0)

keep_species = detected[(detected >= 5) & (var > 0)].index.tolist()
X_raw = X_raw[keep_species]

X_log = np.log1p(X_raw)

X_z_parts = []
for study, idx in crc.groupby("study").groups.items():
    sub = X_log.loc[idx]
    mu = sub.mean(axis=0)
    sd = sub.std(axis=0, ddof=0).replace(0, np.nan)
    z = (sub - mu) / sd
    z = z.fillna(0)
    X_z_parts.append(z)

X = pd.concat(X_z_parts, axis=0).loc[crc.index]
y = crc["stage_binary"].values

dataset_summary = pd.DataFrame([{
    "n_total_crc": len(crc),
    "n_early_stage0_I": int((crc["stage_binary"] == 0).sum()),
    "n_advanced_stageII_IV": int((crc["stage_binary"] == 1).sum()),
    "n_cohorts": crc["study"].nunique(),
    "n_species_input": len(species_cols),
    "n_species_used": len(keep_species),
    "gm_included": "Gemella morbillorum" in keep_species,
}])
dataset_summary.to_csv(res_dir / "01_dataset_summary.csv", index=False)


def cv_feature_importance(model_name, model, X, y):
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    imps = []
    aucs = []

    for train_idx, test_idx in skf.split(X, y):
        X_train = X.iloc[train_idx]
        X_test = X.iloc[test_idx]
        y_train = y[train_idx]
        y_test = y[test_idx]

        model.fit(X_train, y_train)

        if hasattr(model, "predict_proba"):
            score = model.predict_proba(X_test)[:, 1]
        else:
            score = model.decision_function(X_test)

        aucs.append(roc_auc_score(y_test, score))
        imps.append(model.feature_importances_)

    imp = np.mean(np.vstack(imps), axis=0)
    out = pd.DataFrame({
        "species": X.columns,
        "importance": imp,
    }).sort_values("importance", ascending=False).reset_index(drop=True)

    out["rank"] = np.arange(1, len(out) + 1)
    out["method"] = model_name
    out["cv_auc_mean"] = float(np.mean(aucs))
    out["cv_auc_sd"] = float(np.std(aucs))
    return out


extra = cv_feature_importance(
    "Extra Trees",
    ExtraTreesClassifier(
        n_estimators=1000,
        random_state=42,
        class_weight="balanced",
        n_jobs=-1,
    ),
    X,
    y,
)

rf = cv_feature_importance(
    "Random Forest",
    RandomForestClassifier(
        n_estimators=1000,
        random_state=42,
        class_weight="balanced",
        n_jobs=-1,
    ),
    X,
    y,
)

extra_top20 = extra.head(20).copy()
rf_top20 = rf.head(20).copy()

extra.to_csv(res_dir / "01_extra_trees_all_species_cohort_z_log1p.csv", index=False)
rf.to_csv(res_dir / "01_random_forest_all_species_cohort_z_log1p.csv", index=False)

extra_top20.to_csv(
    res_dir / "top20_extra_trees_cohort_z_log1p_stage01_vs_stage234.csv",
    index=False,
)
rf_top20.to_csv(
    res_dir / "top20_random_forest_cohort_z_log1p_stage01_vs_stage234.csv",
    index=False,
)


mw_rows = []
for sp in X.columns:
    early = X_raw.loc[crc["stage_binary"] == 0, sp]
    adv = X_raw.loc[crc["stage_binary"] == 1, sp]

    try:
        stat, p = mannwhitneyu(early, adv, alternative="two-sided")
    except ValueError:
        stat, p = np.nan, 1.0

    mw_rows.append({
        "species": sp,
        "mean_early_raw_abundance": early.mean(),
        "mean_advanced_raw_abundance": adv.mean(),
        "direction": "advanced_high" if adv.mean() > early.mean() else "early_high",
        "mannwhitney_u": stat,
        "p_value": p,
    })

mw = pd.DataFrame(mw_rows)
mw["fdr"] = multipletests(mw["p_value"].fillna(1.0), method="fdr_bh")[1]
mw["minus_log10_fdr"] = -np.log10(mw["fdr"].replace(0, np.nextafter(0, 1)))
mw = mw.sort_values(["fdr", "p_value"]).reset_index(drop=True)
mw["rank"] = np.arange(1, len(mw) + 1)

mw.to_csv(res_dir / "01_mannwhitney_all_species_raw_abundance.csv", index=False)
mw.head(20).to_csv(
    res_dir / "top20_mannwhitney_stage01_vs_stage234.csv",
    index=False,
)

common = sorted(
    set(extra_top20["species"])
    & set(rf_top20["species"])
    & set(mw.head(20)["species"])
)

common_df = pd.DataFrame({"species": common})
common_df["extra_trees_rank"] = common_df["species"].map(extra.set_index("species")["rank"])
common_df["random_forest_rank"] = common_df["species"].map(rf.set_index("species")["rank"])
common_df["mannwhitney_rank"] = common_df["species"].map(mw.set_index("species")["rank"])
common_df.to_csv(res_dir / "01_common_top20_features.csv", index=False)


def plot_top20(df, value_col, title, xlabel, outname):
    d = df.sort_values(value_col, ascending=True)
    plt.figure(figsize=(8, 6))
    plt.barh(d["species"], d[value_col])
    plt.xlabel(xlabel)
    plt.ylabel("")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(fig_dir / f"{outname}.png", dpi=300)
    plt.savefig(fig_dir / f"{outname}.pdf")
    plt.close()


plot_top20(
    extra_top20,
    "importance",
    "Top 20 features: Extra Trees",
    "Feature importance",
    "top20_extra_trees_cohort_z_log1p_stage01_vs_stage234",
)

plot_top20(
    rf_top20,
    "importance",
    "Top 20 features: Random Forest",
    "Feature importance",
    "top20_random_forest_cohort_z_log1p_stage01_vs_stage234",
)

plot_top20(
    mw.head(20),
    "minus_log10_fdr",
    "Top 20 species: Mann-Whitney U test",
    "-log10(FDR)",
    "top20_mannwhitney_stage01_vs_stage234",
)

print("[INFO] Dataset summary")
print(dataset_summary.to_string(index=False))
print("[INFO] Common Top20 features")
print(common_df.to_string(index=False))
print("[INFO] Done: main ML feature selection")
