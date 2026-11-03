# Soil erosion in the Santa River basin upstream of La Balsa (Peru): USLE/RUSLE in Google Earth Engine

Code and derived data for the manuscript *"Formulation uncertainty, paraglacial terrain and mining legacies in water erosion of a deglaciating high-Andean basin: the Santa River upstream of La Balsa (Áncash, Peru)"* (Jesús Aranda, D. R., submitted).

## Contents

| Path | Description |
| --- | --- |
| `notebook/RUSLE_Santa_LaBalsa.ipynb` | Google Colab notebook. Delineates the basin, computes R (4 equations), K, LS (3 formulations), C and P in Google Earth Engine, runs the sensitivity matrix, elevation-belt and deglaciated-terrain analyses, sediment yield and the CHIRPS–gauge comparison. |
| `data/layers/*.png` | Georeferenced RGBA layers exported from Earth Engine (EPSG:4326; bounds in `meta.json`): erosion classes (USLE, RUSLE, RUSLE-MB, RUSLE 2017), hotspots (P90), glaciers 1987–1991 and 2023–2024, deglaciated terrain, R, K, LS, C and drainage network. |
| `data/layers/meta.json` | Layer bounds, P90 threshold, basin polygon, outlet and the 399 mining legacies with RUSLE metrics in a 500 m buffer. |
| `scripts/01_decode_layers.py` | Converts the class and hotspot layers into arrays. |
| `scripts/02_permutation_test.py` | Permutation tests (basin-wide and 2–10 km block displacement of former mining units; 5,000 replicates). |
| `scripts/dem_grid.py` | Resamples the Copernicus GLO-90 DEM to the analysis grid (elevation and slope). |
| `scripts/perm_strat.py` | Permutation test stratified by elevation belt and slope class (5,000 replicates); results in `results/perm_strat_results.json`. |
| `data/external/` | Public monitoring data compiled from ANA (2013, 2016), the Cañón del Pato EIA (1997) and Morera (2010); see `FUENTES.md`. |
| `scripts/03_figures.py` | Reproduces Figures 1–4 of the manuscript. |
| `results/perm_results.json` | Output of the permutation tests. |

## How to run

1. Open the notebook in Google Colab, set `PROJECT_ID` to your Earth Engine project and run all cells (an Earth Engine account is required).
2. For the statistical test and figures (no Earth Engine needed):

```bash
pip install numpy pillow matplotlib
python scripts/01_decode_layers.py
python scripts/02_permutation_test.py
python scripts/03_figures.py
python scripts/dem_grid.py      # downloads Copernicus GLO-90 tiles first (see script)
python scripts/perm_strat.py
```

Random seed: 20261030. Results may differ in the third decimal if the seed is changed.

## Data sources

CHIRPS (Funk et al., 2015), SoilGrids 2.0 (Poggio et al., 2021), NASADEM (NASA JPL, 2020), MERIT Hydro (Yamazaki et al., 2019), HydroBASINS (Lehner and Grill, 2013), Sentinel-2 and Landsat 5 (ESA/USGS), ESA WorldCover v200 (Zanaga et al., 2022) and the national inventory of mining legacies (MINEM, via INGEMMET Geocatmin, accessed 9 October 2026). Only location and classification attributes of the mining legacies are redistributed here.

## License

Code: MIT. Derived data: CC BY 4.0.

## Citation

See `CITATION.cff`.
