# DREADDFUL — DREADDs Fluorescence Ubiety Labeler

*A tiny, transparent helper toolkit around [Cellpose](https://www.cellpose.org/) for quantifying DREADDs expression in the dorsal raphe nucleus (DRN).*  

DREADDFUL is **not** a new segmentation method. It’s a few small, well-documented scripts that:
1) calibrate Cellpose parameters against a subset of manually counted images,  
2) run counting + simple spatial metrics within DRN ROI masks, and  
3) test group differences under cell-count inclusion thresholds.

This repository exists for **methodological transparency** for our manuscript’s (*Behavioral impact of chemogenetic manipulations of 5-HT DRN neurons in transgenic Tph2-iCre rats*) reviewer responses and for anyone wishing to reproduce or adapt our analysis.

> **Ubiety** = “state of being in a particular place” — we use it whimsically to emphasize spatial aspects (coverage / dispersion) of expression inside the DRN mask.

---

## Contents

- `calibrate_cellpose.py` — grid‑search calibration of Cellpose parameters against a small set of manually counted slices.  
- `count_and_spread_cellpose.py` — run Cellpose on all slices, restricted to ROI masks; outputs per‑slice and per‑animal metrics including **n_cells**, **coverage**, and **dispersion**.  
- `group_tests.py` — statistics for group differences under several inclusion thresholds (e.g., ≥50/100/150/200 cells).

---

## Installation (Python)

Tested with **Python 3.10–3.11** and **Cellpose 4.x** (the SAM models). GPU is optional but recommended.

```bash
# create / activate a clean environment (recommended)
python -m venv .venv && source .venv/bin/activate    # (Windows: .venv\\Scripts\\activate)

# install dependencies
pip install "cellpose[all]" scikit-image pandas matplotlib pillow tqdm scikit-learn scipy joblib
```

If you use Conda, create an environment with Python 3.10/3.11 and install the same packages.

---

## Data expectations

### Images & ROI masks
- Each slice image has a **binary ROI mask** that delineates the DRN.  
- By default the scripts look for a mask named `SLICENAME_mask.{tif|tiff|png|jpg}` either
  alongside the image or in `masks/` or `Masks/`.  
- You can also point to masks explicitly via the `mask_path` column in the metadata CSV.

### Metadata CSV
All scripts use a CSV with (at least) these columns:

| column        | required | description |
|---------------|----------|-------------|
| `animal`      | ✓        | animal ID |
| `group`       | ✓        | experimental group label (e.g., `Gi`, `Gq`, `mC`) |
| `slice_path`  | ✓        | absolute or relative path to the image file |
| `mask_path`   |          | optional explicit path to ROI mask (overrides auto‑guess) |
| `manual_count`| (calib)  | manual count used for calibration subset (blank for others) |
| `um_per_px` or `px_per_um` | (calib) | optional pixel size for density readouts |

---

## 1) Calibrate Cellpose parameters

Run a parameter sweep on a **small subset** of manually counted slices and summarize accuracy metrics.

```bash
python calibrate_cellpose.py --config calib_config.json [--gpu]
```

Example `calib_config.json`:

```json
{
  "meta": "path/to/metadata.csv",
  "out_dir": "calib_out",
  "pretrained": "cpsam",
  "threshold": 150,
  "subset": { "slices": ["A1_slice01.tif", "A2_slice09.tif"] },
  "grid": {
    "diameter": [0],
    "cellprob_threshold": [-1.0, -0.5, 0.0],
    "flow_threshold": [0.2, 0.4],
    "min_size": [0, 50, 100],
    "max_size_fraction": [0.0, 0.01],
    "niter": ["auto", 100],
    "tile_overlap": [0.0, 0.1]
  }
}
```

Outputs (inside a timestamped folder):
- `per_slice.csv` — counts & QC per slice for every parameter setting.  
- `per_setting.csv` — summary metrics per setting (e.g., **MAPE**, **balanced accuracy**, **MCC**, **F1**, plus QC medians).  
- `best.json` — top settings by several criteria.

> Our calibration goal was a **≥75% agreement** between auto and manual decisions at the study’s inclusion threshold. (See paper for details.)

---

## 2) Count + spatial metrics on all slices

```bash
python count_and_spread_cellpose.py \
  --meta metadata.csv \
  --out results_dir \
  --config config_cellpose.json \
  --gpu
```

`config_cellpose.json` (global defaults with optional per‑group overrides):

```json
{
  "default": {
    "pretrained_model": "cpsam",
    "diameter": 0,
    "cellprob_threshold": -0.5,
    "flow_threshold": 0.2,
    "normalize": true,

    "min_area": 0,
    "max_area": 0,
    "eccentricity_max": 1.0,
    "solidity_min": 0.0,

    "min_size": 0,
    "max_size_fraction": 0.0,
    "niter": "auto",
    "tile_overlap": 0.0
  },
  "Gi": { "cellprob_threshold": -1.0 }
}
```

Notes:
- `min_size`, `max_size_fraction`, `niter`, and `tile_overlap` are fully supported and recorded in outputs.  
- `max_size_fraction` converts to a per‑slice `max_area = fraction × ROI_area`.  
- The script is defensive to Cellpose API changes: it will drop unsupported kwargs and warn rather than fail.

Key outputs (in a timestamped directory):
- `slice_metrics.csv` — per‑slice **n_cells**, **coverage**, **dispersion** (convex‑hull of centroids / ROI area), plus provenance of key parameters.  
- `centroids.csv` — all detected cell centroids (x,y).  
- `animal_metrics.csv` — slice‑aggregated metrics: `n_cells_total`, `coverage_mean`, `dispersion_mean`.  
- `QC/*.png` — quick overlays of ROI + detected centroids.  
- `labels/*.tif` — label images (post‑filtered within ROI).

---

## 3) Group differences under inclusion thresholds

```bash
python group_tests.py \
  --csv results_dir/DREADDFUL_output_YYMMDD_HHMM/animal_metrics.csv \
  --outdir stats_out \
  --id-col animal \
  --group-col group \
  --count-col n_cells_total \
  --thresholds 50 100 150 200
```

What it does:
- At each threshold T (e.g., **T≥100 cells**), keep animals meeting the count criterion.  
- For **n_cells_total**, **coverage_mean**, **dispersion_mean**:
  - Check **normality** (Shapiro–Wilk, n≥3 per group) and **homoscedasticity** (Levene).  
  - If both pass → **one‑way ANOVA** + Holm‑adjusted t‑tests (Hedges’ g).  
  - Else → **Kruskal–Wallis** + Holm‑adjusted Mann–Whitney U (Cliff’s delta).  
- Saves per‑threshold **descriptives**, **omnibus**, **pairwise**, and **boxplots**.

---

## How this fits the manuscript

- ROI masks were drawn in **Fiji** (ImageJ); 10% of slices were **manually counted** (mCherry+ somas).  
- The calibration minimized **mean absolute/percent error** and targeted **≥75% agreement** with manual decisions at the inclusion threshold.  
- For inclusion we required **≥100** total mCherry+ cells per animal (chosen with reference to prior literature; see manuscript for citations such as *Jones 2022*, *McCosh 2024*, *Radzicki 2024*).  
- Group differences were assessed for **total cell count**, **coverage**, and **dispersion** with the tests described above.

---

## Tips & troubleshooting

- **Masks not found**: ensure each `slice_path` has a matching `*_mask.ext` or set `mask_path` in the CSV.  
- **“not enough valid predictions” (rare)**: try adjusting `cellprob_threshold` upward (stricter) or `flow_threshold`; verify the channel (red).  
- **Large images / memory**: try `tile_overlap > 0` (tiling is auto‑enabled) or downsample images beforehand.  
- **Different Cellpose versions**: the scripts attempt a broad `model.eval(...)` and automatically retry with a conservative subset of kwargs if needed; you’ll see a warning rather than a crash.

---

## Citing & attribution

Please cite:

- **Cellpose** — the core segmentation method.  
  > Pachitariu, M., Rariden, M., & Stringer, C. (2025). Cellpose-SAM: superhuman generalization for cellular segmentation. *bioRxiv.*

- **Fiji/ImageJ** — for ROI mask creation.
  > Schindelin, J., Arganda-Carreras, I., Frise, E., Kaynig, V., Longair, M., Pietzsch, T., … Cardona, A. (2012). Fiji: an open-source platform for biological-image analysis. *Nature Methods*, 9(7), 676–682. doi:10.1038/nmeth.2019

- **This repository (DREADDFUL)** — as a helper toolkit.
  > DREADDFUL (DREADDs Fluorescence Ubiety Labeler): helper scripts for Cellpose‑based quantification of DRN expression. Version 1.0.0. URL: https://github.com/nmccloskey/dreaddful

---

## Acknowledgements

DREADDFUL wraps **Cellpose** and **Fiji**; all credit for segmentation belongs to those projects.

---

## Disclaimer

Built for the DRN chemogenetics project "Behavioral impact of chemogenetic manipulations of 5-HT DRN neurons in transgenic Tph2-iCre rats" by Nicholas S. McCloskey, Chen Li, & Lynn G. Kirby (in review). This is research code provided “as is.” It reflects the exact steps used for the study’s analyses and may need adaptation for other datasets.
