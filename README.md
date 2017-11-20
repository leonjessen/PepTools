<!-- README.md is generated from README.Rmd. Please edit that file -->
PepTools
========

The aim of this package is to supply a set of tools to be used within the fiels of immunoinformatics.

Getting started
---------------

Using Hadley Wickham's brilliant `devtools` package, we can easily install `PepTools` like so:

``` r
install.packages("devtools")
devtools::install_github("leonjessen/PepTools")
```

Once the package has been installed, we can simply load it like so:

``` r
library("PepTools")
#> Loading required package: tidyverse
#> ── Attaching packages ──────────────────────────────────────────────────────────── tidyverse 1.2.1 ──
#> ✔ ggplot2 2.2.1     ✔ purrr   0.2.4
#> ✔ tibble  1.3.4     ✔ dplyr   0.7.4
#> ✔ tidyr   0.7.2     ✔ stringr 1.2.0
#> ✔ readr   1.1.1     ✔ forcats 0.2.0
#> ── Conflicts ─────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::filter() masks stats::filter()
#> ✖ dplyr::lag()    masks stats::lag()
#> Loading required package: ggseqlogo
```

Examples
--------

`PepTools` comes with an example data set of 5,000 9-mer peptides, which have been predicted by [NetMHCpan 4.0 Server](http://www.cbs.dtu.dk/services/NetMHCpan-4.0/) to be strong binders to HLA-A\*02:01.

We can view the first 10 peptides like so

``` r
PEPTIDES %>% head(10)
#>  [1] "YMNSMQEML" "FIYRHMFCV" "VLFKFDMFI" "KLLDRFPVA" "RVLDDFTKL"
#>  [6] "ILKTDILLL" "GLQASIYRI" "AMHAWVFYL" "SMDFGFFEV" "ALPQFNPVV"
```

We can use the `res_counts` function to count the residues:

``` r
PEPTIDES %>% res_counts
#>     A   R   N   D  C   Q   E   G   H   I    L   K    M   F   P   S   T   W
#> 1 324 266 149   8 16 134  24 331 191 322  177 349  282 707   5 365 189 168
#> 2  40   0   0   0  0 264   0   0   0 566 2193   0 1481  14   0  18 106   0
#> 3 286  23 228 409 54 142  80 158 168 366  331  41  501 505 146 191 113 533
#> 4 283 126 317 710 72 219 704 410 223  99   95 214   87 130 479 273 197 128
#> 5 250 151 264 352 53 215 307 388 306 307  200 210  145 306 282 201 180 284
#> 6 226  74 283 266 56 230 199 213 168 515  366 107  286 266 404 260 338 121
#> 7 225  62 219 141 73 209 221  45 342 365  259  58  337 424 443 129 238 420
#> 8 299 241 282 104 46 319 364 233 318 127  238 240  264 312 189 277 379 161
#> 9 444   1   0   0 62   0   1   1   0 872 1148   0  455  64   6  23 127   2
#>     Y    V
#> 1 761  232
#> 2   0  318
#> 3 531  194
#> 4 143   91
#> 5 320  279
#> 6 182  440
#> 7 394  396
#> 8 444  163
#> 9   2 1792
```

Using the `ggseqlogo` package, we can visualise the strong binder peptides in a sequence logo:

``` r
PEPTIDES %>% ggseqlogo()
```

![](README-peps_logo-1.png) ![](README-peps_logo-1.png)
