# Spectral Principal Axis System (SPAS) and tuning of tensor-valued encoding for microscopic anisotropy and time-dependent diffusion

_Code for generating tensor-valued diffusion encoding waveforms, performing simulations, and processing experimental MRI data._
## Overview

This repository contains code supporting the development and evaluation of tensor-valued diffusion MRI encoding strategies to disentangle microscopic anisotropy (µA) from confounding time-dependent diffusion (TDD) effects, published in
[Lasič et al. Spectral Principal Axis System (SPAS) and tuning of tensor-valued encoding for microscopic anisotropy and time-dependent diffusion in the rat brain, Imaging Neuroscience.](https://doi.org/10.3389/fphy.2022.812115)

It implements methods for generating spherical tensor encoding (STE) gradient waveforms and analyzing their spectral properties to define the Spectral Principal Axis System (SPAS). From the same STE, it derives two complementary encoding protocols. The first protocol uses tuned linear tensor encoding (tuned LTE) to minimize TDD bias for unbiased µA estimation. The second protocol uses SPAS-LTE projections, which provide maximal spread of TDD sensitivity inherent in the STE. By acquiring signals along the SPAS-LTE projections and applying geometric averaging of the acquired signals (geoSPAS), this protocol enables both µA estimation and additional contrast sensitive to TDD. The repository provides tools for preparing these encoding waveforms for MRI acquisition and simulations. 
### Citation

If you use this code or methods in your work, please cite:

Lasič, S., Just N., Nilsson, M., Szczepankiewicz, F., Budde, M., & Lundell, H. (2025) Spectral Principal Axis System (SPAS) and tuning of tensor-valued encoding for microscopic anisotropy and time-dependent diffusion in the rat brain. *Imaging Neuroscience*, 9, 812115. https://doi.org/10.3389/fphy.2022.812115

### Requirements

This repository uses the following external toolboxes:
- [Numerical optimization of gradient waveforms (NOW) for tensor-valued dMRI](https://github.com/jsjol/NOW)
- [Multidimensional diffusion MRI repository](https://github.com/markus-nilsson/md-dmri) 
- [Marchenko Pastur PCA denoising repository](https://github.com/sunenj/MP-PCA-Denoising) 
### Code organization

The repository is organized into folders that reflect both functional roles and the general processing flow. Each folder contains executable scripts and supporting functions grouped in subfolders. No centralized script coordinates the workflow. Users are expected to run scripts manually across folders.
### Functional organization

1. **[Waveform generation and preparation](#waveform-generation-and-preparation)** (`make_Bruker_SPAS/`)  
  - Generates initial STE gradient waveforms using numerical optimization (NOW toolbox).
  - Converts SPAS-LTE and tuned LTE waveforms into Bruker-compatible format.
2. **[Spectral analysis, tuning, and simulations](#spectral-analysis-tuning-and-simulations)** (`analyze_simulate_SPAS/`)  
  - Performs spectral analysis of STE waveforms to define the Spectral Principal Axis System (SPAS) and extract SPAS-LTE and tuned LTE projections.
  - Simulates diffusion weighted signals for different waveforms and substrates.
3. **[MRI data analysis](#mri-data-analysis)** (`MRI_Bruker_data_analysis/`)  
  - Processes and analyzes MRI data acquired with the prepared waveforms (Bruker format).
  - Includes scripts for data conversion, denoising, fitting, map generation, visualization, and ROI analysis.
  
### Workflow

The repository supports two complementary workflows: experimental and simulation.

    Generate initial STE waveforms 
          ▼
    Extract SPAS-LTE and tuned LTE projections
          |
          ├──→ Simulate diffusion signals
          ▼
    Convert SPAS-LTE and tuned LTE to Bruker format
          ▼
    Acquire MRI data on scanner
          ▼
    MRI_Bruker_data_analysis/
          ▼
    Process and analyze MRI data

## Documentation

### Waveform generation and preparation
(folder: `make_Bruker_SPAS`)

- `make_wfm1_NOW_optimize.m`  
Generates multiple repetitions of optimized STE gradient waveforms using the [Numerical Optimization of Waveforms (NOW) toolbox](https://github.com/jsjol/NOW).

- `make_wfm2_NOW_interpolate.m`  
Interpolates optimized gradient waveforms, saves the interpolated results, and plots comparisons with the original waveforms.

- `make_wfm3_copy_with_polarity_and_gap.m`  
Copies an interpolated gradient waveform, inserts a gap for the 180° RF pulse, and appends either a polarity-flipped or non-flipped copy.

- `make_wfm3_extra_SE_PGSE.m`  
Generates a PGSE gradient waveform with a 180° RF pulse gap and appends a polarity-flipped copy. 

- `make_wfm4_SPAS2Bruker.m`  
Converts interpolated SPAS-LTE and tuned LTE waveforms into Bruker-compatible gradient waveform files. Generates STE, SPAS-LTE, and LTE outputs, ensures output waveforms match target b-values, and splits waveforms into pre/post blocks. The script expects input waveforms with correct timing and balancing.

### Spectral analysis, tuning, and simulations
(folder: `analyze_simulate_SPAS/`)

- `analyze1_SPAS_tuning_plot.m`  
Main script for analyzing and plotting b-tensor waveforms in multiple steps. Analyze and plot b-tensor waveforms; extract and plot the Spectral Principal Axis System (SPAS) derived as the eigensystem of the low-frequency filtered b-tensor; extract and plot the SPAS-LTE projections; extract and plot the tuned LTE projections for a specified restriction size. Make waveforms `*_SPAS*.mat`, `*_tunedLTE*.mat`, and power spectra `*_info.mat` files.  
  _Functional steps:_
    - step 1: Perform FFT on source b-tensor waveforms  
    - step 2: Show g(t), q(t), spectra, cumulative power trace, b-tensor; extract SPAS and tuned LTE  
    - step 3: Perform FFT on SPAS and tuned LTE waveforms
    - step 4: Display SPAS, tuned waveforms, and power spectra  

  _Spectral anisotropy is visualized using RGB color coding, where the color reflects the relative contributions of each waveform projection within three frequency bands, defined by splitting the cumulative power spectrum at 1/3 and 2/3 of the total encoding power (b-value)._

- `analyze2_tuning_contours_fminsearch.m`  
Generates figures showing: (1) tuning landscapes (proximity to optimal tuning) for varying substrate size and (2) angle vs size between the tuned projections for a selected reference size and other sizes. Uses `fminsearch`; saves outputs to the `fig` folder.

- `analyze3_SPAS_moments_plot.m`  
A helper script to color-code spectral variance (Vw) or centroid frequency and plot the SPAS axis (computed from low-pass filtered b-tensor) alongside spectral moment axes.

- `analyze3_power_moments.m`  
A helper script to check spectral moments such as spectral variance (Vw) or centroid frequency.

- `simulate1_GPA_DvsR_substrates.m`  
Perform frequency domain Gaussian phase approximation (GPA) calculations to simulate apparent diffusion coefficients (ADCs) for various waveforms and substrates (with varying compartment size and shape). Powder averaging is archived by rotating eigen spectra. Requires `*_info.mat` files in `waveform_dir` and combines them into a `wfms_info.mat` file saved in the corresponding `model_*` folder. Generates a `DvsR_*.mat` file containing mean diffusivity (MD) and diffusion variance (VD) as a function of substrate radius for each waveform.

- `simulate2_plot_GPA_DvsR_substrates.m`  
Plot results from GPA simulations stored in `DvsR_*.mat` files generated by `simulate1_GPA_DvsR_substrates.m`. Generate figures showing: (1) signals versus b-value for each waveform and substrate size, (2) signal differences versus size: between geoSPAS and STE, and the maximum difference among the SPAS waveforms, and (3) MD and VD versus size for each waveform.

### MRI data analysis
(folder: `MRI_Bruker_data_analysis/`)

- `step0_setup_code_path.m`  
Initializes MATLAB paths for the MRI Bruker data analysis pipeline. Adds paths to the [Multidimensional diffusion MRI](https://github.com/markus-nilsson/md-dmri) and [MC-PCA denoising](https://github.com/sunenj/MP-PCA-Denoising) repositories.

- `step0_setup_data_path.m`  
Creates a `data_path_struct` structure specifying root data paths, experiment folder names, and ROI file names for Bruker MRI data analysis. Supports defining multiple datasets by adding repeated blocks. Valid paths must be specified in `default_root_data_path`, `exp_folder_name`, and `roi_names` before running data analysis.

- `step1_read_FWF_Bruker.m`  
Reads Bruker datasets, extracts waveform information, and generates NIfTI files and XPS structures (see [Multidimensional diffusion MRI repository](https://github.com/markus-nilsson/md-dmri)) with additional fields for SPAS dataset analysis. Uses the `read_paravision_data` function by Tim Klasen (2013).  

  Requires paths and experiment folders defined in `setup_data_path.m`. Saves output into folders with `_out` suffix unless otherwise specified, including: (1) `info.txt` with scan details, (2) extracted NIfTI files and `_xps.mat` files in `extracted_nii` folders, and (3) waveform information (`*.mat` and `*.png`) in `extracted_wfm` folders. The default setting splits the original data into separate folders for each waveform.

- `step2_merge_Bruker.m`  
Merges NIfTI files and corresponding XPS structures from multiple subfolders into a single merged NIfTI file (`merged_*.nii.gz`) and XPS (`merged_*_xps.mat`). Files can be selected by waveform name or acquisition day; otherwise all valid NIfTI files are merged. Outputs merged files in the experiment folder with prefix `merged_`. Requires input structure from `setup_data_path.m`. Merged NIfTI files can be inspected using `mgui.m` from the [Multidimensional diffusion MRI repository](https://github.com/markus-nilsson/md-dmri).

- `step3_denoise_coreg_Bruker.m`  
Applies [MP-PCA denoising](https://github.com/sunenj/MP-PCA-Denoising) and extrapolation-based motion and eddy-current correction using functions from the [Multidimensional diffusion MRI repository](https://github.com/markus-nilsson/md-dmri). Takes merged NIfTI files (from `step2_merge_Bruker.m`) and outputs denoised and co-registered images. Supports co-registering each dataset independently or merging datasets into a single file for joint co-registration, using Marchenko-Pastur denoising and elastix-based affine registration. See:  
[Veraart et al., NeuroImage (2016) 142, p.394–406](https://doi.org/10.1016/j.neuroimage.2016.08.016),  
[Does et al., Magn Reson Med. 2019; 81:3503–3514](https://doi.org/10.1002/mrm.27658),  
[Nilsson et al., PLoS One. 2015; 10(11):1–22](https://doi.org/10.1371/journal.pone.0141825).

- `step4_processing_Bruker.m`  
Main script for processing Bruker MRI data after merging. Provides a configurable framework for multi-step diffusion MRI analysis, supporting data preprocessing, signal fitting, and parameter mapping, including microscopic anisotropy (µA) and time-dependent diffusion (TDD) contrasts based on raw signal subtraction. 
  
  Input selection is flexible, allowing use of various preprocessed NIfTI files based on specified flags (e.g., `mc` — eddy/motion corrected, `s` — smooth, `pa` — powder average, `g` — geometric average). 

  Outputs include processed NIfTI files and parameter maps, such as:
  - µFA maps (`_muFA_`), 
  - DTI parameter maps (`_dti_`), and 
  - subtraction maps based on raw signal differences (e.g., SPAS1–SPAS3, geoSPAS–STE) for microscopic anisotropy (`_raw_muA_`) and time-dependent diffusion (`_raw_TDD_`) contrasts.

- `step5_make_OP_maps.m`  
Generates order parameter (OP) maps from FA and µFA maps. Takes as input e.g. `'pa_muFA'` and `'dti_fa_geoSPAS'`, outputs NIfTI maps labeled with `_OP_`.

- `step6_subtract_different_maps.m`  
Subtracts different µFA parameter maps to highlight the bias due to time-dependent diffusion (TDD). The µFA maps can be computed in `step4_processing_Bruker.m`, for example using STE and geoSPAS-LTE or STE and SPAS3-LTE, by defining the `opt.SPAS_muFA.LTE_wfm_name`.

- `step7_show_SPAS_ds_maps.m`  
Generates tiled slice images from normalized signal differences and/or log-signal differences at maximum b-value. Maps labeled `_raw_TDD_` and `_raw_muA_` are saved in the `fig` folder. Optionally overlays ROI contours (colored outlines) to highlight selected regions.

- `step7_show_correlation_maps.m`  
Generates joint contrast maps combining subtraction maps µA and TDD (labeled with `_raw_TDD_` and `_raw_muA_`) into a single color-coded visualization. Provides an integrative view of spatial relationships between cell size (TDD) and anisotropy (µA), facilitating identification of tissue regions with distinct microstructural properties. Outputs include tiled slices (saved in the `fig` folder), scatter plots, and optionally segmentation bins (ROIs).

- `step7_show_maps.m`  
Generates tiled slice images (saved in the `fig` folder) from maps produced in previous steps. Can also include subtraction maps from subsequent scans within the same session, generated by the _supplementary script_ `step6_subtract_subsequent_maps.m`.

### Supplementary scripts
(folder: `MRI_Bruker_data_analysis/`) 

These scripts are not part of the main analysis pipeline described in the manuscript.  
They provide additional tools for optional analyses, monitoring, or reporting — useful for exploring changes across subjects or scanning sessions.  The `step#` naming reflects their numbering as organized alongside the main analysis pipeline to suggest the typical order in which they would logically follow and complement the main scripts.

- `step5_subtract_md_SPAS.m`  
Computes difference maps of mean diffusivity (MD) between SPAS1 and SPAS3. Requires input maps labeled `_dti_md_SPAS1` and `_dti_md_SPAS3`. Outputs a difference map labeled `dti_dif_md`, which can optionally be included in ROI statistics (`step6_ROI_stats.m`).

- `step5_ROI_histograms.m`  
Generates overlap histograms and ROI statistics from selected NIfTI maps (e.g., DTI parameter maps, raw signal difference maps) for each ROI and experiment in `select_subfolders`. Saves statistics (mean, standard deviation, confidence intervals) into `stats_all.mat` and histogram figures into the `fig/` subfolder. Useful for comparing ROI data within the same session. 

- `step6_subtract_subsequent_maps.m`  
Subtracts maps from subsequent scans within the same session. Can be used to monitor in vivo changes.

- `step6_ROI_stats.m`  
Performs ROI-based statistics across different experiments for selected NIfTI maps (e.g., DTI parameter maps, µA and TDD subtraction maps, MD difference maps). Computes mean, standard deviation, and confidence intervals for each ROI, map, and experiment. Saves results to `stat.mat`.

- `step7_ROI_stats_make_table.m`  
Creates an Excel table from ROI data saved in `stat.mat` by `step6_ROI_stats.m`. Supports exporting all data or filtered selections by map, ROI, and confidence interval.

- `step8_show_montage.m`  
Assembles montage images by combining saved figure files (e.g., `.png`) from `fig` subfolders across experiments, maps, and figure types (normal maps, subtraction maps, histograms).  Useful for generating compact visual summaries for reporting.

### Function group overview

This section summarizes the main function groups in the repository, listing subfolders as relative paths under each root folder.

#### Functions in `/analyze_simulate_SPAS/functions/`

- `analyze_plot_spectra_SPAS_tuning/`  
  Contains functions for analyzing, plotting, and tuning SPAS spectra, focusing on spectral analysis, visualization, and extraction of the Spectral Principal Axis System (SPAS).

- `Dw_GPA/`  
  Contains functions for frequency-domain Gaussian Phase Approximation (GPA) calculations of ADCs based on tensor-valued diffusion encoding waveforms applied to restricted diffusion in spherical, planar, or cylindrical geometries.

- `SQ/`  
  Contains functions for visualizing tensors using superquadrics.

- `UniformDistSphereRepulsionTri_subset/`  
  Contains precomputed `.mat` files with electrostatic repulsion directional schemes for powder averaging and visualization.

#### Functions in `/MRI_Bruker_data_analysis/functions/`

- `analysis/`  
  Functions supporting diffusion MRI analysis, including:
  - Fitting gamma diffusion distribution and estimating µFA.
  - Computing DTI metrics (MD, FA).
  - Generating parameter maps.
  - Visualizing and analyzing signal decays.

- `color_maps/`  
  Functions for generating perceptually uniform and standardized colormaps, including viridis, plasma, magma, inferno, cividis, and ColorBrewer palettes.

- `image_processing/`  
  Functions for averaging, masking, and smoothing diffusion MRI data, including powder averaging, geometric averaging, repetition averaging, mask creation, and slice-wise Gaussian smoothing.

- `make_maps/`  
  Functions for generating and visualizing maps, including signal difference maps, parameter maps (ADC, FA, OP, µFA), subtraction maps, and correlation maps, with support for ROI analysis.

- `read_Bruker/`  
  Functions for reading and converting Bruker Paravision MRI data and metadata, including `2dseq` file import and preparation of the experimental parameter structure (XPS) for compatibility with the Multidimensional Diffusion MRI repository.

- `utility/`  
  Utility functions supporting MRI data workflows, streamlining tasks such as file handling, data preparation, tensor computations, and visualization.

- `help_scripts/`  
  Contains minor additional helper scripts.
