# Replication code for paper entitled "Strengthened China's food security through air quality improvements".

This directory has three folds:

* /data, which stores the processed and original source data (not shown here because it is empty).
* /figures, which stores the plotted figures from R scripts, but in their specific name.
* /figures_paper, which stores the figures with indexed name.
* /script, which contains the source R code to replicate the results, including processing, preparing, generating results. 

More details in script fold:

* The prefix `load` is for loading packages, functions, and formulas.
* The prefixes `cal_` and `make` demonstrate that they are used for pre-processing source data, for example, how to tidy gridded data to tabular data.
The prefixes `results_` are for generating and plotting the results in Main text and Supporting Information.
* `publish_figs.R` is for rename the figures to figures_paper.
* `tidy_shp.R` is for getting smaller shp files for plotting.
* `tidy_all.R` is for bind the `cal_` data to one tidied dataset.

More details in data folds:

The data in data/. are outputs from `results_` scripts.
The data in data/outputs are outputs from `cal_` scripts.
The data in data/inputs are source data acquired from different sources. Note that due to different data policy, the source data should be download by oneself.

Note that the data fold is empty here due to its large volume. You can contact the authors for acquiring them, if you want.

### R packages

```r
sessioninfo::session_info()
```

```
─ Session info ─────────────────────────────────────────────────────
 setting  value
 version  R version 4.2.1 (2022-06-23)
 os       macOS Ventura 13.1
 system   aarch64, darwin20
 ui       RStudio
 language (EN)
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       Asia/Shanghai
 date     2023-01-11
 rstudio  2022.12.0+353 Elsbeth Geranium (desktop)
 pandoc   2.19.2 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/ (via rmarkdown)

─ Packages ─────────────────────────────────────────────────────────
 package        * version    date (UTC) lib source
 assertthat       0.2.1      2019-03-21 [1] CRAN (R 4.2.0)
 backports        1.4.1      2021-12-13 [1] CRAN (R 4.2.0)
 bayestestR       0.13.0     2022-09-18 [1] CRAN (R 4.2.0)
 biscale        * 1.0.0      2022-05-27 [1] CRAN (R 4.2.0)
 broom          * 1.0.2      2022-12-15 [1] CRAN (R 4.2.0)
 cellranger       1.1.0      2016-07-27 [1] CRAN (R 4.2.0)
 checkmate        2.1.0      2022-04-21 [1] CRAN (R 4.2.0)
 class            7.3-20     2022-01-16 [1] CRAN (R 4.2.1)
 classInt         0.4-8      2022-09-29 [1] CRAN (R 4.2.0)
 cli              3.5.0      2022-12-20 [1] CRAN (R 4.2.0)
 codetools        0.2-18     2020-11-04 [1] CRAN (R 4.2.1)
 collapse       * 1.8.9      2022-10-07 [1] CRAN (R 4.2.0)
 colorspace       2.0-3      2022-02-21 [1] CRAN (R 4.2.0)
 correlation    * 0.8.3      2022-10-09 [1] CRAN (R 4.2.0)
 cowplot        * 1.1.1      2020-12-30 [1] CRAN (R 4.2.0)
 crayon           1.5.2      2022-09-29 [1] CRAN (R 4.2.0)
 data.table     * 1.14.6     2022-11-16 [1] CRAN (R 4.2.0)
 datawizard       0.6.5      2022-12-14 [1] CRAN (R 4.2.0)
 DBI              1.1.3      2022-06-18 [1] CRAN (R 4.2.0)
 dbplyr           2.2.1      2022-06-27 [1] CRAN (R 4.2.0)
 dials          * 1.1.0      2022-11-04 [1] CRAN (R 4.2.1)
 DiceDesign       1.9        2021-02-13 [1] CRAN (R 4.2.0)
 digest           0.6.31     2022-12-11 [1] CRAN (R 4.2.0)
 distributional   0.3.1      2022-09-02 [1] CRAN (R 4.2.0)
 dplyr          * 1.0.10     2022-09-01 [1] CRAN (R 4.2.0)
 dreamerr         1.2.3      2020-12-05 [1] CRAN (R 4.2.0)
 dtplyr         * 1.2.2      2022-08-20 [1] CRAN (R 4.2.0)
 e1071            1.7-12     2022-10-24 [1] CRAN (R 4.2.0)
 ellipsis         0.3.2      2021-04-29 [1] CRAN (R 4.2.0)
 evaluate         0.19       2022-12-13 [1] CRAN (R 4.2.0)
 exactextractr  * 0.9.1      2022-11-16 [1] CRAN (R 4.2.0)
 extrafont        0.18       2022-04-12 [1] CRAN (R 4.2.0)
 extrafontdb      1.0        2012-06-11 [1] CRAN (R 4.2.0)
 fansi            1.0.3      2022-03-24 [1] CRAN (R 4.2.0)
 farver           2.1.1      2022-07-06 [1] CRAN (R 4.2.0)
 fastmap          1.1.0      2021-01-25 [1] CRAN (R 4.2.0)
 fixest         * 0.11.0     2022-10-19 [1] CRAN (R 4.2.0)
 forcats        * 0.5.2      2022-08-19 [1] CRAN (R 4.2.0)
 foreach          1.5.2      2022-02-02 [1] CRAN (R 4.2.0)
 Formula          1.2-4      2020-10-16 [1] CRAN (R 4.2.0)
 fs               1.5.2      2021-12-08 [1] CRAN (R 4.2.0)
 furrr          * 0.3.1      2022-08-15 [1] CRAN (R 4.2.0)
 future         * 1.30.0     2022-12-16 [1] CRAN (R 4.2.0)
 future.apply     1.10.0     2022-11-05 [1] CRAN (R 4.2.1)
 gargle           1.2.1      2022-09-08 [1] CRAN (R 4.2.0)
 gdtools          0.2.4      2022-02-14 [1] CRAN (R 4.2.0)
 generics         0.1.3      2022-07-05 [1] CRAN (R 4.2.0)
 ggdensity      * 0.1.1      2022-10-24 [1] CRAN (R 4.2.0)
 ggdist         * 3.2.0      2022-07-19 [1] CRAN (R 4.2.0)
 ggh4x          * 0.2.3      2022-11-09 [1] CRAN (R 4.2.1)
 gghighlight    * 0.4.0      2022-10-16 [1] CRAN (R 4.2.0)
 ggnewscale       0.4.8      2022-10-06 [1] CRAN (R 4.2.0)
 ggplot2        * 3.4.0      2022-11-04 [1] CRAN (R 4.2.1)
 ggpmisc        * 0.5.2      2022-12-17 [1] CRAN (R 4.2.0)
 ggpointdensity * 0.1.0      2019-08-28 [1] CRAN (R 4.2.0)
 ggpp           * 0.5.0      2022-12-05 [1] CRAN (R 4.2.0)
 ggrepel        * 0.9.2      2022-11-06 [1] CRAN (R 4.2.0)
 ggside         * 0.2.2      2022-12-04 [1] CRAN (R 4.2.0)
 ggstar         * 1.0.4      2022-11-08 [1] CRAN (R 4.2.1)
 ggupset        * 0.3.0      2020-05-05 [1] CRAN (R 4.2.0)
 globals          0.16.2     2022-11-21 [1] CRAN (R 4.2.0)
 glue             1.6.2      2022-02-24 [1] CRAN (R 4.2.0)
 googledrive      2.0.0      2021-07-08 [1] CRAN (R 4.2.0)
 googlesheets4    1.0.1      2022-08-13 [1] CRAN (R 4.2.0)
 gower            1.0.1      2022-12-22 [1] CRAN (R 4.2.0)
 GPfit            1.0-8      2019-02-08 [1] CRAN (R 4.2.0)
 gridExtra        2.3        2017-09-09 [1] CRAN (R 4.2.0)
 gtable           0.3.1      2022-09-01 [1] CRAN (R 4.2.0)
 hardhat          1.2.0      2022-06-30 [1] CRAN (R 4.2.0)
 haven            2.5.1      2022-08-22 [1] CRAN (R 4.2.0)
 hms              1.1.2      2022-08-19 [1] CRAN (R 4.2.0)
 hrbrthemes     * 0.8.0      2020-03-06 [1] CRAN (R 4.2.0)
 htmltools        0.5.4      2022-12-07 [1] CRAN (R 4.2.0)
 httr             1.4.4      2022-08-17 [1] CRAN (R 4.2.0)
 infer          * 1.0.4      2022-12-02 [1] CRAN (R 4.2.0)
 insight          0.18.8     2022-11-24 [1] CRAN (R 4.2.0)
 ipred            0.9-13     2022-06-02 [1] CRAN (R 4.2.0)
 iterators        1.0.14     2022-02-05 [1] CRAN (R 4.2.0)
 janitor          2.1.0      2021-01-05 [1] CRAN (R 4.2.0)
 jsonlite         1.8.4      2022-12-06 [1] CRAN (R 4.2.0)
 KernSmooth       2.23-20    2021-05-03 [1] CRAN (R 4.2.1)
 knitr            1.41       2022-11-18 [1] CRAN (R 4.2.1)
 latex2exp      * 0.9.6      2022-11-28 [1] CRAN (R 4.2.0)
 lattice          0.20-45    2021-09-22 [1] CRAN (R 4.2.1)
 lava             1.7.1      2023-01-06 [1] CRAN (R 4.2.1)
 lazyeval         0.2.2      2019-03-15 [1] CRAN (R 4.2.0)
 lhs              1.1.6      2022-12-17 [1] CRAN (R 4.2.0)
 lifecycle        1.0.3      2022-10-07 [1] CRAN (R 4.2.0)
 listenv          0.9.0      2022-12-16 [1] CRAN (R 4.2.0)
 lubridate      * 1.9.0      2022-11-06 [1] CRAN (R 4.2.0)
 magrittr       * 2.0.3      2022-03-30 [1] CRAN (R 4.2.0)
 MASS             7.3-58.1   2022-08-03 [1] CRAN (R 4.2.0)
 Matrix         * 1.5-3      2022-11-11 [1] CRAN (R 4.2.0)
 MatrixModels     0.5-1      2022-09-11 [1] CRAN (R 4.2.0)
 matrixStats    * 0.63.0     2022-11-18 [1] CRAN (R 4.2.0)
 MetBrewer        0.2.0      2022-03-21 [1] CRAN (R 4.2.0)
 mgcv           * 1.8-41     2022-10-21 [1] CRAN (R 4.2.0)
 modeldata      * 1.0.1      2022-09-06 [1] CRAN (R 4.2.0)
 modelr           0.1.10     2022-11-11 [1] CRAN (R 4.2.1)
 munsell          0.5.0      2018-06-12 [1] CRAN (R 4.2.0)
 mvtnorm          1.1-3      2021-10-08 [1] CRAN (R 4.2.0)
 ncdf4          * 1.20       2022-12-03 [1] CRAN (R 4.2.0)
 nlme           * 3.1-161    2022-12-15 [1] CRAN (R 4.2.0)
 nnet             7.3-18     2022-09-28 [1] CRAN (R 4.2.0)
 numDeriv         2016.8-1.1 2019-06-06 [1] CRAN (R 4.2.0)
 pammtools      * 0.5.8      2022-01-09 [1] CRAN (R 4.2.0)
 parallelly       1.33.0     2022-12-14 [1] CRAN (R 4.2.0)
 parsnip        * 1.0.3      2022-11-11 [1] CRAN (R 4.2.1)
 patchwork      * 1.1.2      2022-08-19 [1] CRAN (R 4.2.0)
 pec              2022.05.04 2022-05-04 [1] CRAN (R 4.2.0)
 pillar           1.8.1      2022-08-19 [1] CRAN (R 4.2.0)
 pkgconfig        2.0.3      2019-09-22 [1] CRAN (R 4.2.0)
 prodlim          2019.11.13 2019-11-17 [1] CRAN (R 4.2.0)
 proxy            0.4-27     2022-06-09 [1] CRAN (R 4.2.0)
 purrr          * 1.0.0      2022-12-20 [1] CRAN (R 4.2.0)
 purrrgress     * 0.0.1      2022-09-25 [1] Github (tylergrantsmith/purrrgress@b9c4daa)
 qs             * 0.25.4     2022-08-09 [1] CRAN (R 4.2.0)
 quantreg         5.94       2022-07-20 [1] CRAN (R 4.2.0)
 R6               2.5.1      2021-08-19 [1] CRAN (R 4.2.0)
 RApiSerialize    0.1.2      2022-08-25 [1] CRAN (R 4.2.0)
 raster         * 3.6-11     2022-11-28 [1] CRAN (R 4.2.0)
 Rcpp             1.0.9      2022-07-08 [1] CRAN (R 4.2.0)
 RcppParallel     5.1.5      2022-01-05 [1] CRAN (R 4.2.0)
 readr          * 2.1.3      2022-10-01 [1] CRAN (R 4.2.0)
 readxl         * 1.4.1      2022-08-17 [1] CRAN (R 4.2.0)
 recipes        * 1.0.3      2022-11-09 [1] CRAN (R 4.2.1)
 reprex           2.0.2      2022-08-17 [1] CRAN (R 4.2.0)
 rlang            1.0.6      2022-09-24 [1] CRAN (R 4.2.0)
 rmarkdown        2.19       2022-12-15 [1] CRAN (R 4.2.0)
 rpart            4.1.19     2022-10-21 [1] CRAN (R 4.2.0)
 rsample        * 1.1.1      2022-12-07 [1] CRAN (R 4.2.0)
 rstudioapi       0.14       2022-08-22 [1] CRAN (R 4.2.0)
 Rttf2pt1         1.3.11     2022-10-08 [1] CRAN (R 4.2.0)
 rvest            1.0.3      2022-08-19 [1] CRAN (R 4.2.0)
 sandwich         3.0-2      2022-06-15 [1] CRAN (R 4.2.0)
 scales         * 1.2.1      2022-08-20 [1] CRAN (R 4.2.0)
 sessioninfo      1.2.2      2021-12-06 [1] CRAN (R 4.2.0)
 sf             * 1.0-9      2022-11-08 [1] CRAN (R 4.2.0)
 showtext         0.9-5      2022-02-09 [1] CRAN (R 4.2.0)
 showtextdb       3.0        2020-06-04 [1] CRAN (R 4.2.0)
 snakecase        0.11.0     2019-05-25 [1] CRAN (R 4.2.0)
 sp             * 1.5-1      2022-11-07 [1] CRAN (R 4.2.0)
 SparseM          1.81       2021-02-18 [1] CRAN (R 4.2.0)
 stringfish       0.15.7     2022-04-13 [1] CRAN (R 4.2.0)
 stringi          1.7.8      2022-07-11 [1] CRAN (R 4.2.0)
 stringr        * 1.5.0      2022-12-02 [1] CRAN (R 4.2.0)
 survival         3.4-0      2022-08-09 [1] CRAN (R 4.2.0)
 sysfonts         0.8.8      2022-03-13 [1] CRAN (R 4.2.0)
 systemfonts      1.0.4      2022-02-11 [1] CRAN (R 4.2.0)
 terra          * 1.6-47     2022-12-02 [1] CRAN (R 4.2.0)
 tibble         * 3.1.8      2022-07-22 [1] CRAN (R 4.2.0)
 tidymodels     * 1.0.0      2022-07-13 [1] CRAN (R 4.2.0)
 tidyr          * 1.2.1      2022-09-08 [1] CRAN (R 4.2.0)
 tidyselect       1.2.0      2022-10-10 [1] CRAN (R 4.2.0)
 tidyverse      * 1.3.2      2022-07-18 [1] CRAN (R 4.2.0)
 timechange     * 0.1.1      2022-11-04 [1] CRAN (R 4.2.0)
 timeDate         4022.108   2023-01-07 [1] CRAN (R 4.2.1)
 timereg          2.0.4      2022-11-09 [1] CRAN (R 4.2.0)
 tune           * 1.0.1      2022-10-09 [1] CRAN (R 4.2.0)
 tzdb             0.3.0      2022-03-28 [1] CRAN (R 4.2.0)
 units            0.8-1      2022-12-10 [1] CRAN (R 4.2.0)
 utf8             1.2.2      2021-07-24 [1] CRAN (R 4.2.0)
 vctrs            0.5.1      2022-11-16 [1] CRAN (R 4.2.0)
 viridis        * 0.6.2      2021-10-13 [1] CRAN (R 4.2.0)
 viridisLite    * 0.4.1      2022-08-22 [1] CRAN (R 4.2.0)
 withr            2.5.0      2022-03-03 [1] CRAN (R 4.2.0)
 workflows      * 1.1.2      2022-11-16 [1] CRAN (R 4.2.1)
 workflowsets   * 1.0.0      2022-07-12 [1] CRAN (R 4.2.0)
 WrensBookshelf   0.1.0      2022-08-15 [1] CRAN (R 4.2.0)
 xfun             0.36       2022-12-21 [1] CRAN (R 4.2.0)
 xml2             1.3.3      2021-11-30 [1] CRAN (R 4.2.0)
 yardstick      * 1.1.0      2022-09-07 [1] CRAN (R 4.2.0)
 zoo              1.8-11     2022-09-17 [1] CRAN (R 4.2.0)

 [1] /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library

────────────────────────────────────────────────────────────────────
```
