# China_SIF_air_pollution_Code

This is the replication code for paper entitled 'Strengthened China's food security through air quality improvements'.

This directory has three folds:

1. data, which stores the processed and original source data (not shown here because it is empty).
2. figures, which stores the plotted figures from R scripts, but in their specific name.
3. figures_paper, which stores the figures with indexed name.
4. script, which contains the source R code to replicate the results, including processing, preparing, generating results. 

More details in script fold:

1. The prefix `load` is for loading packages, functions, and formulas.
2. The prefixes `cal_` and `make` demonstrate that they are used for pre-processing source data, for example, how to tidy gridded data to tabular data.
The prefixes `results_` are for generating and plotting the results in Main text and Supporting Information.
4. `publish_figs.R` is for rename the figures to figures_paper.
5. `tidy_shp.R` is for getting smaller shp files for plotting.
6. `tidy_all.R` is for bind the `cal_` data to one data.

More details in data folds:

The data in data/. are outputs from `results_` scripts.
The data in data/outputs are outputs from `cal_` scripts.
The data in data/inputs are source data acquired from different sources. Note that due to different data policy, the source data should be download by oneself.

Note that the data fold does not contains due to its large volume. You can contact the authors for acquiring them, if you want.
