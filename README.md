# Effective Connectivity of the Human Insula

## About

This repository contains scripts and data used to generate matrices and brain maps illustrating results from:

**_Effective connectivity of the human insula measured by cortico-cortical evoked potentials_**  
Pinheiro et al.

The analysis is based on stereo-EEG (sEEG) stimulation data from the F-TRACT database and focuses on effective connectivity between the insula and whole-brain regions.

The latest version of the data will be made available through **EBRAINS** (currently in progress).

---

## Repository Structure

### Folders

- **results/**  
  Contains all generated data and figures. See *Results Folder Structure* below.  

- **parcellation_definitions/**  
  Labels and definitions for each supported parcellation scheme.

- **mne_data/**  
  MNE-based resources used for cortical surface visualization.

- **ENIGMA/**  
  ENIGMA toolbox resources used for subcortical mesh visualization.

### Scripts

- **analyze_atlas.py**  
  Connectivity analysis workflow for a single atlas, including statistical analysis, visualization of brain maps, and insula dominance analysis.

- **compare_atlas.py**  
  Performs statistical comparisons between two atlases and generates directional connectivity brain maps.

- **dicROIs_Insula.py**  
  Parcel indices and labels for insula low, mid, mid-to-high, and high resolution according to Montreal (MNI-insula) and Julich parcellation schemes.

### Environment

- **environment2.yml**  
  Conda environment file containing all required dependencies.

---

## Results Folder Structure

### Data Organization

Each data folder follows a fixed hierarchical structure that encodes the analysis configuration.

Example path:\
\ft-insula\data\results\INSULA\
\Lausanne2018-scale1\nan\MNI-insula\19\efferent_False\
\max_peak_delay_100.0__zth5\
\feature_ampl_zth5__N_with_value__min_value__50

The folder hierarchy is defined as follows (from top to bottom):

1. **Stimulation parcellation**  
   Parcellation scheme used for stimulation (e.g., `Lausanne2018`, `MNI-insula`, `MNI-JulichBrain-3.0`).

2. **Stimulation parcellation resolution**  
   Number of parcels included in the stimulation parcellation.  
   A value of `nan` indicates the raw (native) parcellation resolution.

3. **Recording parcellation**  
   Parcellation scheme used for recording sites.

4. **Recording parcellation resolution**  
   Number of parcels included in the recording parcellation.  
   A value of `nan` indicates the raw resolution.

5. **Connectivity direction**  
   Indicates whether connectivity is analyzed as:
   - `efferent_True`  (insula stimulation → recording whole brain)
   - `efferent_False` (stimulation whole brain → recording insula)

6. **Temporal constraint and z-threshold**  
   Folder encoding the maximum peak delay (in ms) and z-score threshold used to detect responses.

7. **Minimum response requirement**  
Folder encoding the feature type, z-threshold, and minimum required number of responses.

---

### Data Format

Each configuration contains **parcel-by-parcel matrices** for multiple features.

- **Rows** correspond to **stimulation parcels**
- **Columns** correspond to **recording parcels**

#### Amplitude Features
- Total number of responses (`NTotal`)
- Number of significant responses (`N_with_value`)
- Descriptive statistics of first peak amplitude (z-score):
  - mean, median (absolute)
  - min, max
  - standard deviation
  - quantiles (0.25, 0.5, 0.75)

#### Probability
- Number of significant responses (`N_with_value`) divided by the total number of responses (`NTotal`)

#### Peak Delay Features
- Number of significant responses
- Descriptive statistics of first peak latency (ms)

#### Implantation Features
- Number of implantations (`count_unique_str`)

#### Patient Features
- Number of patients (`count_unique_str`)

---

### Available Stimulation × Recording Parcellation Combinations

The repository includes multiple combinations of stimulation and recording parcellation schemes, supporting both **efferent** and **afferent** connectivity analyses.

Efferent Connectivity *(Insula stimulation → whole-brain recording)*:

- **Montreal (MNI-insula; 1, 4, 7, 19 parcels)** × Lausanne2018-scale2
- **Julich (1, 4, 10, 16 parcels)** × Lausanne2018-scale2
- **Montreal / Julich (19 / 16 parcels)** × Lausanne2018-scale1
- **Montreal / Julich (4 / 4 parcels)** × Lausanne2018-scale3

Afferent Connectivity *(Whole-brain stimulation → insula recording)*:

- Lausanne2018-scale2 × **Montreal (1, 4, 7, 19 parcels)**
- Lausanne2018-scale2 × **Julich (1, 4, 10, 16 parcels)**
- Lausanne2018-scale1 × **Montreal / Julich (19 / 16 parcels)**
- Lausanne2018-scale3 × **Montreal / Julich (4 / 4 parcels)**

---

## analyze_atlas.py

### Description

Connectivity analysis workflow for a single atlas, including statistical analysis, visualization of brain maps, and insula dominance analysis.

### Workflow Overview

1. **Data paths and parcel definitions**
   - Identify and define all available data paths
   - Define stimulation and recording parcels (supported parcellations: Lausanne, Montreal, and Julich)
   - Update parcel names according to optional parcel merging (aggregates data across all parcels within a merged group) and hemispheric symmetrization (aggregates data across inter- and intra-hemispheric connections)

2. **Data loading**
   - Discrete measures:
     ```
     probability, NTotal, NValue, patients, implantations
     ```
   - Continuous measures:
     ```
     ampl, peak_delay
     ```
   - Continuous features:
     ```
     mean, std, quantile_0.25, quantile_0.5, quantile_0.75
     ```

3. **Parcel inclusion mask**
   - Filter parcels based on minimum response count and minimum number of patients
   - Compute and visualize brain coverage and confidence intervals of probabilistic connectivity after applying the inclusion mask

4. **Matrices and brain maps**
   - Matrices plotted for each data type
   - Brain maps generated for efferent or afferent connectivity
   - Cortical parcels visualized using FreeSurfer meshes
   - Subcortical parcels:
     - Lausanne: ENIGMA template meshes
     - Other parcellations: 2D subcortical plots

5. **Mesh visualization**
   - Left and right hemispheres (when hemispheric symmetrization is applied, intra- and inter-hemispheric connections are displayed as left and right, respectively)
   - Lateral and medial views
   - Meshes and associated labels saved in mne_data folder as:
        ```
        lh.<parcel_name>.label
        rh.<parcel_name>.label
        ```
     
6. **Insula dominance analysis (only for resolution > 1)**
   - Perform pairwise proportion z-tests between insula connections targeting the same brain parcel
   - Effect size computed using Cohen’s *h* 
   - Benjamini–Hochberg correction
   - Dominance matrix:
     - `2`: highest probability parcel
     - `1`: statistically similar to the highest
     - `0`: remaining insula parcels
   - Generate a bar plot showing the percentage of dominant insula parcels per brain system (occipital, temporal, parietal, frontal, limbic, central) based on their connectivity with the insula. Brain systems only defined for Lausanne parcellation. 

---

### Inputs (analyze_atlas.py)

| Parameter | Description                                                                                                               |
|---------|---------------------------------------------------------------------------------------------------------------------------|
| `resolution_stim` | Stimulation parcellation resolution (e.g., 1)                                                                             |
| `resolution_rec` | Recording parcellation resolution (e.g., np.nan)                                                                          |
| `stim_parc_name` | Stimulation parcellation scheme (e.g., "MNI-insula")                                                                      |
| `rec_parc_name` | Recording parcellation scheme (e.g., "Lausanne2018-scale2")                                                               |
| `plot_parc_name` | Parcellation used for visualization (e.g., "Lausanne2018-scale2")                                                         |
| `merged` | Parcel merging configuration for stimulation and recording parcellations (e.g., [True, False])                            |
| `hemis_symmetrize` | Hemisphere symmetrization (`True` or `False`) and matrix layout (`"H"` for horizontal, `"V"` for vertical) (e.g., [True, "H"]) |
| `ntotal_min` | Minimum required responses (e.g., 50)                                                                                     |
| `time_window` | Peak delay limit (e.g., `["max", 100]` means the peak delay is considered until a maximum of 100 ms)                      |                                                                                         |
| `efferent` | Efferent (`True`) or afferent (`False`) connectivity (e.g., True)                                                         |
| `save` | Save figures (`True or False`)                                                                                            |

### Outputs (analyze_atlas.py)

Saved as `.svg` files:

- **Matrices:** probability, NTotal, patients, implantations, amplitude, peak delay, dominance
- **Brain maps:** one map per insula parcel

---

## compare_atlas.py

### Description

Performs statistical comparisons between probabilistic connectivity of two atlases and generates directional connectivity brain maps.

### Supported Comparisons

- Efferent vs afferent atlases (EffvsAff_)
- Left vs right hemisphere (LeftvsRight_)
- All contacts vs gray-matter-only contacts (AllMatter_vs_Gray_)

### Workflow Overview

1. **Data paths and parcel definitions**
2. **Data loading**
3. **Parcel inclusion mask**
4. **Matrices and brain maps**
5. **Mesh visualization**
6. **Pairwise differences between two atlases (only supports probability data)**
   - Proportion z-tests for each corresponding connection between Atlas 1 and Atlas 2 
   - Directionality is computed by subtracting the probabilities of the corresponding connections between Atlas 1 and Atlas 2 (Atlas 1 − Atlas 2) and normalizing the result to the range −1 to +1 using the maximum probability
   - Benjamini–Hochberg correction 

### Inputs

- `title` – Type of comparison between atlases. Options include: `EffvsAff_`, `LeftvsRight_`, `AllMatter_vs_Gray_`

- Same parameters as `analyze_atlas.py`, specified for **two atlases** (e.g., `[[True, "H"], [True, "V"]]` for hemis_symmetrize input).

### Outputs

- Directionality brain maps (`.svg`)
- Figures saved in the efferent atlas results directory

---

## Citation

If you use this code or data, cite as follows:

Pinheiro et al., *Effective connectivity of the human insula measured by cortico-cortical evoked potentials*


