# Overview

Repository for the analysis of electrophysiological data after spike sorting. Analyses are focused to single cell medial prefrontal cortex (mPFC) activity during an auditory delayed comparison working memory task (Parametric Working Memory, PMW) in rats. Moreover, they are inspired by those done by [Romo ... Brody et al. 1999](https://pubmed.ncbi.nlm.nih.gov/10365959/) in macaques.

## Highlights
* functions for aligning `kilosort` spike sorting output with `Bdata` SQL tables containing behavioral data
* flexible python code for exploratory visualizations of neural activity during stimulus and delay period conditioned on multiple features
* inspired by the [Spykes](https://github.com/jess-breda/spykes) library & rewritten with updated python syntax & seaborn visuals

![data_sdc_20190902_145404_fromSD_N0_delay_overlap_dual_delay_plot](https://user-images.githubusercontent.com/53059059/202951294-ef9dc87a-1808-47f3-9869-0ad10aeb782d.png)

## Steps
- (1) Align spike and behavioral data (`\kilosort_bdata_alignment`)
- (2) Import to python (`io_utils.py`)
- (3) Visualization and analyses (`plotting_utils.py`)

## Usage
```
git clone  https://github.com/Brody-Lab/jbreda_PWM_mPFC_analysis
```
```
conda create -n PWM_ephys python=3.7 pip numpy matplotlib scipy scikit-learn h5py pyqt cython pillow black seaborn jupyter statsmodels pydove
```
```
conda activate PWM_ephys
pip3 install jedi==0.17.2 # for tab complete bug
```

### Tutorial

See the `2021_Tutorial.ipynb` notebook for example usage.


### Assumptions

(1) This repo is built as a continuation of the output from the [jbreda_kilsort](https://github.com/Brody-Lab/jbreda_kilosort) repository that was designed for analyzing wireless tetrode data in `kilosort` and uses the same file structure.

(2) Functions for spike & behavior alignment are written to be run *after* all sorting is done on a per animal basis. But, if needed, this could be easily changed

