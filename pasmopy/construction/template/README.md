# Template for BioMASS model construction

A brief description of each file/folder is below:

| Name                                           | Content                                                                                                  |
| ---------------------------------------------- | -------------------------------------------------------------------------------------------------------- |
| [`name2idx/`](./name2idx/)                     | Names of model parameters and species                                                                    |
| [`set_model.py`](./set_model.py)               | Differential equation, parameters and initial condition                                                  |
| [`observalbe.py`](./observable.py)             | Observables, simulations and experimental data                                                           |
| [`viz.py`](./viz.py)                           | Plotting parameters for customizing figure properties                                                    |
| [`set_search_param.py`](./set_search_param.py) | Model parameters to optimize and search region                                                           |
| [`fitness.py`](./fitness.py)                   | An objective function to be minimized, i.e., the distance between model simulation and experimental data |
| [`reaction_network.py`](./reaction_network.py) | Reaction indices grouped according to biological processes                                               |

## Examples

- [Circadian clock](https://github.com/biomass-dev/biomass/tree/master/biomass/models/circadian_clock)
- [Immediate-early gene response](https://github.com/biomass-dev/biomass/tree/master/biomass/models/Nakakuki_Cell_2010)
- [Insulin signaling](https://github.com/biomass-dev/biomass/tree/master/biomass/models/insulin_signaling)
- [MAPK cascade](https://github.com/biomass-dev/biomass/tree/master/biomass/models/mapk_cascade)
- [NF-κB pathway](https://github.com/biomass-dev/biomass/tree/master/biomass/models/nfkb_pathway)
- [TGF-β/SMAD pathway](https://github.com/biomass-dev/biomass/tree/master/biomass/models/tgfb_smad)