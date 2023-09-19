# Replication code for paper entitled "Strengthening China's food security through air quality improvements".

This directory has three folds:

* /data, which stores the processed and original source data (not shown here because it is empty).
* /figures, which stores the plotted figures from R scripts, but in their specific name.
* /figures_paper, which stores the figures with indexed name.
* /script, which contains the source R code to replicate the results, including processing, preparing, and generating data and figures. 

More details in script fold:

* The prefix `load` is for loading packages, functions, and formulas.
* The prefixes `cal_` and `make` demonstrate that they are used for pre-processing source data, for example, how to tidy gridded data to tabular data.
The prefixes `results_` are for generating and plotting the results in Main text and Supporting Information.
* `publish_figs.R` is for rename the figures in /figures to /figures_paper.
* `tidy_shp.R` is for getting smaller shp files for plotting.
* `tidy_all.R` is for bind the `cal_` data to one tidied dataset.

More details in data fold:

The data in data/. are outputs from `results_` scripts.
The data in data/outputs are outputs from `cal_` scripts.
The data in data/inputs are source data acquired from different sources. Note that due to different data policy, the source data should be download by oneself.

Note that the data fold is empty here due to its large volume. Please contact the authors for acquiring them (email: xliu21@smail.nju.edu.cn), if you want.

### R packages

```r
sessioninfo::session_info()
```

```
─ Session info ────────────────────────────────────────────────────────────
 setting  value
 version  R version 4.2.3 (2023-03-15 ucrt)
 os       Windows 10 x64 (build 19045)
 system   x86_64, mingw32
 ui       RStudio
 language (EN)
 collate  Chinese (Simplified)_China.utf8
 ctype    Chinese (Simplified)_China.utf8
 tz       Asia/Taipei
 date     2023-09-19
 rstudio  2022.12.0+353 Elsbeth Geranium (desktop)
 pandoc   2.19.2 @ D:/Software/RStudio/resources/app/bin/quarto/bin/tools/ (via rmarkdown)

─ Packages ────────────────────────────────────────────────────────────────
 ! package           * version    date (UTC) lib source
   backports           1.4.1      2021-12-13 [1] CRAN (R 4.2.0)
   bayestestR          0.13.1     2023-04-07 [1] CRAN (R 4.2.3)
   biscale           * 1.0.0      2022-05-27 [1] CRAN (R 4.2.3)
   broom             * 1.0.5      2023-06-09 [1] CRAN (R 4.2.3)
   cellranger          1.1.0      2016-07-27 [1] CRAN (R 4.2.3)
   class               7.3-22     2023-05-03 [1] CRAN (R 4.2.3)
   classInt            0.4-10     2023-09-05 [1] CRAN (R 4.2.3)
   cli                 3.6.1      2023-03-23 [1] CRAN (R 4.2.3)
   cluster             2.1.4      2022-08-22 [1] CRAN (R 4.2.3)
   codetools           0.2-19     2023-02-01 [1] CRAN (R 4.2.3)
   collapse          * 1.9.6      2023-05-28 [1] CRAN (R 4.2.3)
   colorspace          2.1-0      2023-01-23 [1] CRAN (R 4.2.3)
   correlation       * 0.8.4      2023-04-06 [1] CRAN (R 4.2.3)
   cowplot           * 1.1.1      2020-12-30 [1] CRAN (R 4.2.3)
   crayon              1.5.2      2022-09-29 [1] CRAN (R 4.2.3)
   crul                1.4.0      2023-05-17 [1] CRAN (R 4.2.3)
   curl                5.0.2      2023-08-14 [1] CRAN (R 4.2.3)
   data.table        * 1.14.8     2023-02-17 [1] CRAN (R 4.2.3)
   datawizard          0.9.0      2023-09-15 [1] CRAN (R 4.2.3)
   DBI                 1.1.3      2022-06-18 [1] CRAN (R 4.2.3)
   deldir              1.0-9      2023-05-17 [1] CRAN (R 4.2.3)
   dials             * 1.2.0      2023-04-03 [1] CRAN (R 4.2.3)
   DiceDesign          1.9        2021-02-13 [1] CRAN (R 4.2.3)
   digest              0.6.33     2023-07-07 [1] CRAN (R 4.2.3)
   distributional      0.3.2      2023-03-22 [1] CRAN (R 4.2.3)
   dplyr             * 1.1.3      2023-09-03 [1] CRAN (R 4.2.3)
   dreamerr            1.3.0      2023-08-23 [1] CRAN (R 4.2.3)
   dtplyr            * 1.3.1      2023-03-22 [1] CRAN (R 4.2.3)
   e1071               1.7-13     2023-02-01 [1] CRAN (R 4.2.3)
   ellipsis            0.3.2      2021-04-29 [1] CRAN (R 4.2.3)
   evaluate            0.21       2023-05-05 [1] CRAN (R 4.2.3)
   exactextractr     * 0.9.1      2022-11-16 [1] CRAN (R 4.2.3)
   extrafont           0.19       2023-01-18 [1] CRAN (R 4.2.2)
   extrafontdb         1.0        2012-06-11 [1] CRAN (R 4.2.0)
   fansi               1.0.4      2023-01-22 [1] CRAN (R 4.2.3)
   farver              2.1.1      2022-07-06 [1] CRAN (R 4.2.3)
   fastmap             1.1.1      2023-02-24 [1] CRAN (R 4.2.3)
   fixest            * 0.11.1     2023-01-10 [1] CRAN (R 4.2.3)
   fontBitstreamVera   0.1.1      2017-02-01 [1] CRAN (R 4.2.0)
   fontLiberation      0.1.0      2016-10-15 [1] CRAN (R 4.2.0)
   fontquiver          0.2.1      2017-02-01 [1] CRAN (R 4.2.3)
   forcats           * 1.0.0      2023-01-29 [1] CRAN (R 4.2.3)
   foreach             1.5.2      2022-02-02 [1] CRAN (R 4.2.3)
   Formula             1.2-5      2023-02-24 [1] CRAN (R 4.2.2)
   furrr             * 0.3.1      2022-08-15 [1] CRAN (R 4.2.3)
   future            * 1.33.0     2023-07-01 [1] CRAN (R 4.2.3)
   future.apply        1.11.0     2023-05-21 [1] CRAN (R 4.2.3)
   gdtools             0.3.3      2023-03-27 [1] CRAN (R 4.2.3)
   generics            0.1.3      2022-07-05 [1] CRAN (R 4.2.3)
   gfonts              0.2.0      2023-01-08 [1] CRAN (R 4.2.3)
   ggdensity         * 1.0.0      2023-02-09 [1] CRAN (R 4.2.3)
   ggdist            * 3.3.0      2023-05-13 [1] CRAN (R 4.2.3)
   ggh4x             * 0.2.6      2023-08-30 [1] CRAN (R 4.2.3)
   gghighlight       * 0.4.0      2022-10-16 [1] CRAN (R 4.2.3)
   ggplot2           * 3.4.3      2023-08-14 [1] CRAN (R 4.2.3)
   ggpmisc           * 0.5.4-1    2023-08-13 [1] CRAN (R 4.2.3)
   ggpointdensity    * 0.1.0      2023-04-02 [1] Github (LKremer/ggpointdensity@02f3ab2)
   ggpp              * 0.5.4      2023-08-12 [1] CRAN (R 4.2.3)
   ggrepel           * 0.9.3      2023-02-03 [1] CRAN (R 4.2.3)
   ggside            * 0.2.2      2022-12-04 [1] CRAN (R 4.2.3)
   ggstar            * 1.0.4      2022-11-08 [1] CRAN (R 4.2.3)
   ggupset           * 0.3.0      2020-05-05 [1] CRAN (R 4.2.3)
   globals             0.16.2     2022-11-21 [1] CRAN (R 4.2.2)
   glue                1.6.2      2022-02-24 [1] CRAN (R 4.2.3)
   gower               1.0.1      2022-12-22 [1] CRAN (R 4.2.2)
   GPfit               1.0-8      2019-02-08 [1] CRAN (R 4.2.3)
   gridExtra           2.3        2017-09-09 [1] CRAN (R 4.2.3)
   gtable              0.3.4      2023-08-21 [1] CRAN (R 4.2.3)
   hardhat             1.3.0      2023-03-30 [1] CRAN (R 4.2.3)
   hexbin              1.28.3     2023-03-21 [1] CRAN (R 4.2.3)
   hms                 1.1.3      2023-03-21 [1] CRAN (R 4.2.3)
   hrbrthemes        * 0.8.0      2020-03-06 [1] CRAN (R 4.2.3)
   htmltools           0.5.6      2023-08-10 [1] CRAN (R 4.2.3)
   httpcode            0.3.0      2020-04-10 [1] CRAN (R 4.2.3)
   httpuv              1.6.11     2023-05-11 [1] CRAN (R 4.2.3)
   infer             * 1.0.5      2023-09-06 [1] CRAN (R 4.2.3)
   insight             0.19.5     2023-09-13 [1] CRAN (R 4.2.3)
   interp              1.1-4      2023-03-31 [1] CRAN (R 4.2.3)
   ipred               0.9-14     2023-03-09 [1] CRAN (R 4.2.3)
   iterators           1.0.14     2022-02-05 [1] CRAN (R 4.2.3)
   jpeg                0.1-10     2022-11-29 [1] CRAN (R 4.2.2)
   jsonlite            1.8.7      2023-06-29 [1] CRAN (R 4.2.3)
   KernSmooth          2.23-22    2023-07-10 [1] CRAN (R 4.2.3)
   knitr               1.44       2023-09-11 [1] CRAN (R 4.2.3)
   later               1.3.1      2023-05-02 [1] CRAN (R 4.2.3)
   latex2exp         * 0.9.6      2022-11-28 [1] CRAN (R 4.2.3)
   lattice             0.21-8     2023-04-05 [1] CRAN (R 4.2.3)
   latticeExtra        0.6-30     2022-07-04 [1] CRAN (R 4.2.3)
   lava                1.7.2.1    2023-02-27 [1] CRAN (R 4.2.3)
   lhs                 1.1.6      2022-12-17 [1] CRAN (R 4.2.3)
   lifecycle           1.0.3      2022-10-07 [1] CRAN (R 4.2.3)
   listenv             0.9.0      2022-12-16 [1] CRAN (R 4.2.3)
   lubridate         * 1.9.2      2023-02-10 [1] CRAN (R 4.2.3)
   magrittr          * 2.0.3      2022-03-30 [1] CRAN (R 4.2.3)
   mapproj             1.2.11     2023-01-12 [1] CRAN (R 4.2.3)
   maps                3.4.1      2022-10-30 [1] CRAN (R 4.2.3)
   MASS                7.3-60     2023-05-04 [1] CRAN (R 4.2.3)
   Matrix            * 1.6-1      2023-08-14 [1] CRAN (R 4.2.3)
   MatrixModels        0.5-2      2023-07-10 [1] CRAN (R 4.2.3)
   matrixStats       * 1.0.0      2023-06-02 [1] CRAN (R 4.2.3)
   mgcv              * 1.9-0      2023-07-11 [1] CRAN (R 4.2.3)
   mime                0.12       2021-09-28 [1] CRAN (R 4.2.0)
   modeldata         * 1.2.0      2023-08-09 [1] CRAN (R 4.2.3)
   munsell             0.5.0      2018-06-12 [1] CRAN (R 4.2.3)
   ncdf4             * 1.21       2023-01-07 [1] CRAN (R 4.2.2)
   nlme              * 3.1-163    2023-08-09 [1] CRAN (R 4.2.3)
   nnet                7.3-19     2023-05-03 [1] CRAN (R 4.2.3)
   numDeriv            2016.8-1.1 2019-06-06 [1] CRAN (R 4.2.0)
   openair           * 2.17-0     2023-05-02 [1] CRAN (R 4.2.3)
   parallelly          1.36.0     2023-05-26 [1] CRAN (R 4.2.3)
   parsnip           * 1.1.1      2023-08-17 [1] CRAN (R 4.2.3)
   patchwork         * 1.1.3      2023-08-14 [1] CRAN (R 4.2.3)
   pillar              1.9.0      2023-03-22 [1] CRAN (R 4.2.3)
   pkgconfig           2.0.3      2019-09-22 [1] CRAN (R 4.2.3)
   png                 0.1-8      2022-11-29 [1] CRAN (R 4.2.2)
   polynom             1.4-1      2022-04-11 [1] CRAN (R 4.2.3)
   prodlim             2023.08.28 2023-08-28 [1] CRAN (R 4.2.3)
   promises            1.2.1      2023-08-10 [1] CRAN (R 4.2.3)
   proxy               0.4-27     2022-06-09 [1] CRAN (R 4.2.3)
   purrr             * 1.0.2      2023-08-10 [1] CRAN (R 4.2.3)
   qs                * 0.25.5     2023-02-22 [1] CRAN (R 4.2.3)
   quantreg            5.97       2023-08-19 [1] CRAN (R 4.2.3)
   R6                  2.5.1      2021-08-19 [1] CRAN (R 4.2.3)
   RApiSerialize       0.1.2      2022-08-25 [1] CRAN (R 4.2.1)
   raster            * 3.6-23     2023-07-04 [1] CRAN (R 4.2.3)
   RColorBrewer        1.1-3      2022-04-03 [1] CRAN (R 4.2.0)
   Rcpp                1.0.11     2023-07-06 [1] CRAN (R 4.2.3)
 D RcppParallel        5.1.7      2023-02-27 [1] CRAN (R 4.2.3)
   readr             * 2.1.4      2023-02-10 [1] CRAN (R 4.2.3)
   readxl            * 1.4.3      2023-07-06 [1] CRAN (R 4.2.3)
   recipes           * 1.0.8      2023-08-25 [1] CRAN (R 4.2.3)
   rlang               1.1.1      2023-04-28 [1] CRAN (R 4.2.3)
   rmarkdown           2.25       2023-09-18 [1] CRAN (R 4.2.3)
   rpart               4.1.19     2022-10-21 [1] CRAN (R 4.2.3)
   rsample           * 1.2.0      2023-08-23 [1] CRAN (R 4.2.3)
   rstudioapi          0.15.0     2023-07-07 [1] CRAN (R 4.2.3)
   Rttf2pt1            1.3.12     2023-01-22 [1] CRAN (R 4.2.2)
   sandwich            3.0-2      2022-06-15 [1] CRAN (R 4.2.3)
   scales            * 1.2.1      2022-08-20 [1] CRAN (R 4.2.3)
   sessioninfo         1.2.2      2021-12-06 [1] CRAN (R 4.2.3)
   sf                * 1.0-14     2023-07-11 [1] CRAN (R 4.2.3)
   shiny               1.7.5      2023-08-12 [1] CRAN (R 4.2.3)
   showtext            0.9-6      2023-05-03 [1] CRAN (R 4.2.3)
   showtextdb          3.0        2020-06-04 [1] CRAN (R 4.2.3)
   slider            * 0.3.0      2022-11-16 [1] CRAN (R 4.2.3)
   sp                * 2.0-0      2023-06-22 [1] CRAN (R 4.2.3)
   SparseM             1.81       2021-02-18 [1] CRAN (R 4.2.0)
   stringfish          0.15.8     2023-05-30 [1] CRAN (R 4.2.3)
   stringi             1.7.12     2023-01-11 [1] CRAN (R 4.2.2)
   stringr           * 1.5.0      2022-12-02 [1] CRAN (R 4.2.3)
   survival            3.5-7      2023-08-14 [1] CRAN (R 4.2.3)
   sysfonts            0.8.8      2022-03-13 [1] CRAN (R 4.2.3)
   systemfonts         1.0.4      2022-02-11 [1] CRAN (R 4.2.3)
   terra             * 1.7-46     2023-09-06 [1] CRAN (R 4.2.3)
   tibble            * 3.2.1      2023-03-20 [1] CRAN (R 4.2.3)
   tidymodels        * 1.1.1      2023-08-24 [1] CRAN (R 4.2.3)
   tidyr             * 1.3.0      2023-01-24 [1] CRAN (R 4.2.3)
   tidyselect          1.2.0      2022-10-10 [1] CRAN (R 4.2.3)
   tidyverse         * 2.0.0      2023-02-22 [1] CRAN (R 4.2.3)
   timechange          0.2.0      2023-01-11 [1] CRAN (R 4.2.3)
   timeDate            4022.108   2023-01-07 [1] CRAN (R 4.2.3)
   tune              * 1.1.2      2023-08-23 [1] CRAN (R 4.2.3)
   tzdb                0.4.0      2023-05-12 [1] CRAN (R 4.2.3)
   units               0.8-3      2023-08-10 [1] CRAN (R 4.2.3)
   utf8                1.2.3      2023-01-31 [1] CRAN (R 4.2.3)
   vctrs               0.6.3      2023-06-14 [1] CRAN (R 4.2.3)
   viridis           * 0.6.4      2023-07-22 [1] CRAN (R 4.2.3)
   viridisLite       * 0.4.2      2023-05-02 [1] CRAN (R 4.2.3)
   warp                0.2.0      2020-10-21 [1] CRAN (R 4.2.3)
   withr               2.5.0      2022-03-03 [1] CRAN (R 4.2.3)
   workflows         * 1.1.3      2023-02-22 [1] CRAN (R 4.2.3)
   workflowsets      * 1.0.1      2023-04-06 [1] CRAN (R 4.2.3)
   xfun                0.40       2023-08-09 [1] CRAN (R 4.2.3)
   xtable              1.8-4      2019-04-21 [1] CRAN (R 4.2.3)
   yardstick         * 1.2.0      2023-04-21 [1] CRAN (R 4.2.3)
   zoo                 1.8-12     2023-04-13 [1] CRAN (R 4.2.3)

────────────────────────────────────────────────────────────────────
```
