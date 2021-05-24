## Parameter estimation

Using experimental datasets obtained from breast cancer cell lines to identify model parameters.

### Requirements

| Language                            | Package                                              | Desctiption                                                                                                           |
| ----------------------------------- | ---------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------- |
| [Julia 1.5+](https://julialang.org) | [`BioMASS.jl`](https://github.com/himoto/BioMASS.jl) | This module provides a Julia interface to the [biomass](https://github.com/okadalabipr/biomass) parameter estimation. |

### Run optimization

```bash
$ mkdir errout
$ sh optimize_parallel.sh
```

### How to track optimization process

The temporary result will be saved in [`erbb_network_jl/logs/n.log`](https://github.com/pasmopy/breast_cancer/tree/master/training/erbb_network_jl/logs) after each iteration.

```bash
$ tail erbb_network_jl/logs/1.log
```
