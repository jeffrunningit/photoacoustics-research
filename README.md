# Laser-Induced Photoacoustics in Aluminium Gratings

## Overview

This project explores laser-induced photoacoustics in a 130 nm aluminium grating, focusing on the effects of surface plasmon polaritons (SPPs) on the detection of acoustic waves. The study combines experimental data collection, computational analysis, and visualization to understand laser-induced strain wave propagation in the grating.

## Directory Structure

- **Data/**  
  Contains experimental data files, processing scripts, and analysis routines.  
  - **AFM/**: AFM measurements of the grating surface (height and phase data).  
  - **Resonance angle calculation data/**:  
    - *`Al resonance angle calculation.py`*: Calculates the resonance angle for aluminium.  
    - *`plot Al res wl vs angle.py`*: Plots the resonance wavelength versus incidence angle.  
  - **grating 750nm measurement/**:  
    Raw measurement data organized by date (e.g., 20240212, 20240213, etc.). Each date folder contains:  
    - Raw data files (.obt, .dat, .txt)  
    - *`data_processing_script.py`*: Processes the raw data (averaging, delta calculations, etc.) and generates outputs.  
  - **grating spectrum/**:  
    Contains spectroscopic data and scripts such as *`spectrometer plot.py`* and *`spectrum_plot_script.py`* for plotting reflection spectra.  
  - **grating whitelight measurement/**:  
    - *`dataprocessing_jeffrey.py`*: Processes whitelight measurement data.  
    - *`whitelight pumpprobe 2dplot.py`*: Generates 2D pump–probe plots from the data.  
  - **Additional scripts**:  
    - *`peak fluence.py`*: Calculates or processes the laser peak fluence.  
    - *`data into npzfile.py`*: Converts raw data into NPZ format.  
    - *`grating 750nm probe measurement plots.py`*: Visualizes probe measurement data.  
  - **750nm -1 to +1 data** and **750nm -3 to +3 data**:  
    Include scripts (*plot_script.py*, *data_processing_script.py*) for analyzing data taken at different angular configurations.

- **Figures/**  
  Contains generated figures and plots (PNG, PDF, SVG) that visualize experimental results, such as AFM cross-sectional images, transient grating measurements, and spectroscopic data. Some plotting scripts (e.g., *crosssection_plot.py*, *transient_grating_plot.py*) are also found here.

- **Thesis/**  
  Contains the final thesis document and presentation that detail the experimental work and analysis.


## Code Files Description

- **Al resonance angle calculation.py**: Computes the resonance angle based on experimental/theoretical parameters for the aluminium grating.
- **plot Al res wl vs angle.py**: Plots the resonance wavelength versus incidence angle.
- **data_processing_script.py**: Found in various measurement directories, these scripts process raw data (e.g., averaging, delta analysis) to prepare it for visualization.
- **spectrometer plot.py** & **spectrum_plot_script.py**: Visualize spectroscopic data from the grating measurements.
- **dataprocessing_jeffrey.py** & **whitelight pumpprobe 2dplot.py**: Handle the processing and plotting of whitelight measurement data.
- **peak fluence.py**: Calculates or processes laser peak fluence.
- **data into npzfile.py**: Converts raw experimental data into NPZ files for streamlined analysis.
- **grating 750nm probe measurement plots.py**: Processes and plots probe measurement data.
- **plot_script.py**, **crosssection_plot.py**, **transient_grating_plot.py**: Scripts for generating various plots and visualizations from the experimental data.

## Requirements

- **Python 3.x**  
- Python packages:  
  - NumPy  
  - SciPy  
  - Matplotlib  
  - Pandas  
  - (Other packages may be required; please check individual scripts for additional dependencies.)

## Usage

1. **Data Processing**:  
   - Navigate to the specific data subdirectory (e.g., `Data/grating 750nm measurement/20240212`).  
   - Run the corresponding processing script:  
     ```bash
     python data_processing_script.py
     ```
   - Processed data (e.g., averaged results, delta values) and output figures will be generated in the same folder.

2. **Plot Generation**:  
   - To generate plots (such as resonance angle plots or spectroscopic graphs), run the associated plotting script, for example:  
     ```bash
     python plot_Al_res_wl_vs_angle.py
     ```
   - Figures will be saved to the designated output directory.

3. **Data Conversion**:  
   - If needed, run the utility script to convert raw data into NPZ format:  
     ```bash
     python data_into_npzfile.py
     ```

4. **Whitelight Measurements**:  
   - Process and plot whitelight data using:  
     ```bash
     python dataprocessing_jeffrey.py
     python whitelight_pumpprobe_2dplot.py
     ```

## Notes

- Ensure that all required dependencies are installed before running the scripts.
- The experimental data are organized by date and measurement type. Verify that you are working in the correct directory for the data set of interest.
- Some scripts are optimized for specific datasets and may require minor modifications for use with other data.

## Acknowledgements

This project is part of an in-depth study on laser-induced photoacoustics and plasmonic effects in aluminium gratings. The detailed experimental work and analysis are documented in the Thesis folder.