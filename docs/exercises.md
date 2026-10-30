---
title: Plotting exercises
description: Fetch GW and ZTF data and make example plots for the Hands-on school
permalink: /exercises/
---

Hands-on notebooks that **fetch** official public data and **plot** example figures for the ZTF–GW / AGN-channel project. Run each with `uv` from the repo root after `uv sync`.

```bash
uv run jupyter notebook notebooks/exercises/
# or one file:
uv run jupyter notebook notebooks/exercises/01_fetch_plot_gw.ipynb
```

## Exercise 1 — GW data

**Notebook:** [`notebooks/exercises/01_fetch_plot_gw.ipynb`](https://github.com/AstroSuro/EAGW2026/blob/main/notebooks/exercises/01_fetch_plot_gw.ipynb)

1. Fetch open **H1/L1 strain** around GW190521 with GWpy / GWOSC and plot a band-passed time series.
2. Fetch the official **PE sky-map FITS** with `scripts/fetch_data.py` and plot a Mollweide map with the 90% contour.

```bash
uv run python scripts/fetch_data.py skymaps --events GW190521
```

## Exercise 2 — ZTF flare data

**Notebook:** [`notebooks/exercises/02_fetch_plot_ztf.ipynb`](https://github.com/AstroSuro/EAGW2026/blob/main/notebooks/exercises/02_fetch_plot_ztf.ipynb)

1. Fetch the IRSA light curve for **ZTF19abanrhr**.
2. Plot magnitude vs MJD by filter and mark the GW190521 time.

```bash
uv run python scripts/fetch_data.py ztf --object ZTF19abanrhr
```

## Exercise 3 — ZTF position on the GW sky map

**Notebook:** [`notebooks/exercises/03_overlay_ztf_on_skymap.ipynb`](https://github.com/AstroSuro/EAGW2026/blob/main/notebooks/exercises/03_overlay_ztf_on_skymap.ipynb)

Overlay the AGN / ZTF19abanrhr coordinates on the GW190521 PE sky map and report the credible level that encloses that pixel — a first spatial look at the flare–GW association.

## Prep checklist

| Step | Command / notebook |
|---|---|
| Environment | `uv sync` |
| First light | `uv run jupyter notebook notebooks/first_light.ipynb` |
| Optional bulk fetch | `uv run python scripts/fetch_data.py all --events GW190521` |
| Exercises 1–3 | `notebooks/exercises/01_*.ipynb` … `03_*.ipynb` |

Finish first light and at least Exercises 1–2 by **13 November**.
