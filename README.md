# Oceanic Spatio-temporal Patterns of Vertical Velocities in the Cape Basin and Agulhas Current Retroflection From Two Years of SWOT Altimetry

Authors: S. Coadou-Chaventon, L. Siegelman<sup>\*</sup>, E. Carli<sup>\*</sup>, S. Swart, R. Schubert, S. Speich

This is a repository for the data analyses and figures of *Oceanic Spatio-temporal Patterns of Vertical Velocities in the Cape Basin and Agulhas Current Retroflection From Two Years of SWOT Altimetry* (in review for GRL). 

### Abstract
<p align="justify">
The Surface Water and Ocean Topography (SWOT) satellite provides a groundbreaking view of sea surface height across a 120-km-wide swath (20-km nadir gap), opening new opportunities to reconstruct vertical velocities ($w$). Vertical motions play a central role in redistributing properties, influencing climate and ecosystems. Using the effective surface quasigeostrophic framework in the eddy-rich Agulhas Current Retroflection region, we reconstruct $w$ down to 1,000 m. We assess its realism against the 1/60° (~2 km) INALT60 ocean model, finding a spatial correlation of 0.6. SWOT-derived $w$ reveal pronounced high-frequency variability, with events lasting several days to one week during which $w$ variance more than doubles. The Agulhas Retroflection emerges as a hotspot of enhanced $w$ (~300 m day<sup>-1</sup>), driven by strong horizontal strain. This work provides the first quantitative reconstruction of $w$ from SWOT in the region and highlights the mission’s potential to quantify vertical exchanges linking the ocean surface and interior.
</p>

### Workflow

The notebooks are organized as follows:

1. Derivation of daily root mean square (RMS) profiles of $w$ from INALT60 over the entire simulation period. Used to optimize the parameters for the effective Surface Quasi Geostrophy (eSQG), N<sup>2</sup> and c - [01_RMS_w_INALT60.ipynb](https://github.com/SolangeCoadou/CoadouChaventon-2026-GRL/blob/main/01_RMS_w_INALT60.ipynb)
2. Tuning of N<sup>2</sup> and c for each month on the ten days of INALT60 most representative of the monthly-mean vertical profile of the vertical velocity RMS - [02_Optimization_N2_and_c.ipynb](https://github.com/SolangeCoadou/CoadouChaventon-2026-GRL/blob/main/02_Optimization_N2_and_c.ipynb)
3. Statistics on the derived N<sup>2</sup> and c. Shows no clear seasonal cycle in the parameters value but strong variability within the monthly N<sup>2</sup> distributions - [03_Timeseries_N2_and_c.ipynb](https://github.com/SolangeCoadou/CoadouChaventon-2026-GRL/blob/main/03_Timeseries_N2_and_c.ipynb)
4. Processing INALT60 fields used in Figure 1 - [04_Fields_INALT60_extraction.ipynb](https://github.com/SolangeCoadou/CoadouChaventon-2026-GRL/blob/main/04_Fields_INALT60_extraction.ipynb)
5. Plotting Figure 1. Shows the region of interest and some INALT60 filtered fields. It also highlights the moderate correlation between the vertical velocity variance and the strain magnitude - [05_Plot_INALT60_fields.ipynb](https://github.com/SolangeCoadou/CoadouChaventon-2026-GRL/blob/main/05_Plot_INALT60_fields.ipynb)
6. Evaluating the eSQG performance in deriving relative vorticity and vertical velocity fields from the sea surface height (SSH) using INALT60. Reveals that the eSQG framework robustly reconstruct w associated with mesoscale features - [06_Correlation_eSQG_INALT60.ipynb](https://github.com/SolangeCoadou/CoadouChaventon-2026-GRL/blob/main/06_Correlation_eSQG_INALT60.ipynb)
7. Reconstruction of $w$ over the SWOT fast sampling phase. Shows strong regional variability high-frequency variability, with events lasting several days to one week during which $w$ variance more than doubles - [07_Fast_sampling_reconstruction.ipynb](https://github.com/SolangeCoadou/CoadouChaventon-2026-GRL/blob/main/07_Fast_sampling_reconstruction.ipynb)
8. Processing INALT60 fields to derive the properties shown in Figure 4 - [08_SSH_w_and_OW_INALT60.ipynb](https://github.com/SolangeCoadou/CoadouChaventon-2026-GRL/blob/main/08_SSH_w_and_OW_INALT60.ipynb)
9. Computing statistics from INALT60 vertical velocity field shown in Figure 4  - [09_INALT60_statistics.ipynb](https://github.com/SolangeCoadou/CoadouChaventon-2026-GRL/blob/main/09_INALT60_statistics.ipynb)
10. Reconstruction of $w$ over the SWOT science phase. Reveals the Agulhas Retroflection as a hotspot of intense $w$, resulting from the strong horizontal strain - [10_Science_phase_reconstruction.ipynb](https://github.com/SolangeCoadou/CoadouChaventon-2026-GRL/blob/main/10_Science_phase_reconstruction.ipynb)

All the functions are available in a separate file (functions.py) and loaded at the beginning of each of the notebooks.

Some of the code used for solving the effective Surface Quasi Geostrophy equations has been developped by Elisa Carli, Lia Siegelman and Patrice Klein [(Carli et al., 2025)](https://zenodo.org/records/15088480). The spectral computations rely on the work of [Vivant et al., 2025](https://www.nature.com/articles/s43247-025-02002-z), available [here](https://zenodo.org/records/13923050).   
</p> 

### Data sources

* [SWOT altimetry](https://www.aviso.altimetry.fr/en/data/products/sea-surface-height-products/global/swot-l3-ocean-products.html)
* [DUACS altimetry](https://data.marine.copernicus.eu/product/SEALEVEL_GLO_PHY_L4_MY_008_047/description)
* [INALT60 numerical model simulation](https://hdl.handle.net/20.500.12085/58605cbf-99da-4a9f-a070-98e346230583), Schubert & Schwarzkopf, 2026


