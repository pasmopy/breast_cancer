---
layout: default
title: Requirements
nav_order: 2
---

# Requirements

## Manual installation of package requirements

The following packages are required for transcriptomic data integration, parameter estimation, patient-specific simulations, and result visualization.

### General:

- Python >= 3.7
- Julia >= 1.5
- R >= 4.0

### Python:

```
$ pip install -r requirements.txt
```

> - [pasmopy==0.1.0](https://github.com/pasmopy/pasmopy)
> - [biomass>=0.5.2,<0.6](https://github.com/biomass-dev/biomass)
> - [matplotlib==3.3.4](https://matplotlib.org)
> - [numpy==1.19.2](https://numpy.org)
> - [pandas==1.2.4](https://pandas.pydata.org)
> - [seaborn==0.11.2](https://seaborn.pydata.org)
> - [scipy==1.6.0](https://scipy.org)

### Julia:

```shell
$ julia
julia> ] add BioMASS @0.5.0
```

> - [BioMASS.jl==0.5.0](https://github.com/biomass-dev/BioMASS.jl)

### R:

```shell
$ R
> source("install_requirements.R")
```

> - [TCGAbiolinks](https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html)
> - [sva](https://bioconductor.org/packages/release/bioc/html/sva.html)
> - [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)
> - [ComplexHeatmap](https://www.bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)
> - [viridisLite](https://github.com/sjmgarnier/viridisLite)
> - [dplyr](https://dplyr.tidyverse.org)
> - [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
> - [tibble](https://tibble.tidyverse.org)
> - [data.table](https://github.com/Rdatatable/data.table)
> - [stringr](https://stringr.tidyverse.org)
