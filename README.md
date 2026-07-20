# netCoin: Interactive Analytic Networks, Galleries and Plots for Coincidences and Regressions

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/netCoin)](https://CRAN.R-project.org/package=netCoin)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/netCoin)](https://CRAN.R-project.org/package=netCoin)
[![License: GPL (>= 2)](https://img.shields.io/badge/license-GPL%20(%3E%3D%202)-blue.svg)](https://www.gnu.org/licenses/gpl-2.0)
<!-- badges: end -->

**netCoin** combines the data analysis capabilities of R with the 
interactive visualization libraries of JavaScript to create networks of coincidences, 
co-occurrences, correlations and regressions that can be explored directly 
in a browser, embedded in HTML files or integrated into Shiny applications. 
It also generates networked HTML galleries and other plots of coincidences.

## Mission

This project's aim is to integrate traditional statistical techniques with
automatic learning and social network analysis tools for the purpose of
obtaining visual and interactive displays of big data. The interdisciplinary
team involved has the following objectives:

1. Efficiently combine different statistical techniques by integrating them
   under the study of the coincidence of people, objects, events or
   characteristics in a multiple series of scenarios.
2. Design open-source software that, under the premise of network coincidence
   analysis, generates different types of interactive graphics that enable an
   exploratory and confirmatory analysis to be made of vast quantities of
   information.
3. Apply all the above to the creation and handling of large databases in such
   diverse fields as the following:
   - survey data combined with administrative data;
   - the analysis of networks created by Twitter users and those reproduced
     through their messages;
   - the abstracts of scientific output in different disciplines over long
     periods of time through the generation of semantic maps;
   - the creation of a huge database of leading figures in the fields of
     philosophy, science, social sciences and the arts, which also contains
     their major works.

## Installation

Install the released version from CRAN:

```r
install.packages("netCoin")
```

Or the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("Modesto-Escobar/netCoin")
```

## Quick start

```r
library(netCoin)

# Classic coincidence network from a binary incidence matrix
data(finches)
net <- netCoin(finches)
plot(net)                      # opens an interactive HTML visualisation

# Coincidence analysis with statistical significance
co <- coin(finches)
summary(co)
plot(co)
```

## Main features

| Area | Key functions |
|---|---|
| General coincidence networks | `netCoin()`, `addNetCoin()`, `allNet()` |
| Classic coincidence analysis | `coin()`, `coocur()`, `propCoin()` |
| Correlation networks | `netCorr()`, `d_netCorr()` |
| Survey-oriented networks | `surCoin()`, `surScat()` |
| Regression-based networks | `logCoin()`, `glmCoin()` |
| Path and cobweb networks | `pathCoin()`, `cobCoin()` |
| Galleries and multi-panel views | `gallery()`, `netGallery()`, `netExhibit()`, `multiPages()` |
| Layout and node handling | `asNodes()`, `layoutCircle()`, `layoutGrid()` |
| Export to external formats | `savePajek()`, `saveGhml()` |
| Shiny integration | `shinyCoin()` |

Interactive visualisations are powered by D3.js through the companion package
[`rD3plot`](https://CRAN.R-project.org/package=rD3plot).

## Documentation

- Package website: <https://modesto-escobar.github.io/netCoin/>
- Vignettes (after installation):

  ```r
  vignette("netCoin")     # introduction and main use cases
  vignette("surCoin")     # survey-oriented networks
  vignette("galleries")   # building galleries and exhibits
  ```

## Citation

To cite `netCoin` in publications, run:

```r
citation("netCoin")
```

## Authors

- **Modesto Escobar** (Universidad de Salamanca) — creator & maintainer
  [[ORCID]](https://orcid.org/0000-0003-2072-6071)
- David Barrios (Universidad de Salamanca)
- Carlos Prieto (Universidad de Salamanca)
  [[ORCID]](https://orcid.org/0000-0003-2064-4842)
- Luis Martínez-Uribe (Universidad de Salamanca)
  [[ORCID]](https://orcid.org/0000-0002-7795-3972)
- Pablo Cabrera-Álvarez (University of Essex)
  [[ORCID]](https://orcid.org/0000-0001-8105-5908)
- Cristina Calvo-López (Universidad Nacional de Educación a Distancia)
  [[ORCID]](https://orcid.org/0000-0001-5039-1263)

## Acknowledgments

This work has been supported by grants **CSO2013-49278-EXP**,
**PGC2018-093755-B100**, **PDC2022-133355-100** and **PID2023-147358NB-100**,
funded by MICIU/AEI/10.13039/501100011033 and by the European Union
NextGenerationEU/PRTR programme.

## Bug reports and contributions

Please file issues and pull requests at:

<https://github.com/Modesto-Escobar/netCoin-2.x/issues>

## License

GPL-2 | GPL-3
