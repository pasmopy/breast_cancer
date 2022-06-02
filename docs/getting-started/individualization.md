---
layout: default
title: Individualization of the mechanistic model
parent: Getting started
nav_order: 3
---

# Individualization of the mechanistic model

## Parameter estimation

Here, we use phospho-protein time-course datasets to train kinetic constants and weighting factors.

📒 **Note:** In this study, we used `BioMASS.jl`, but in most cases you can use `pasmopy.optimize()` function for parameter estimation.

1. Build a mechanistic model to identify model parameters

   ```python
   import os

   from pasmopy import Text2Model


   Text2Model(os.path.join("models", "erbb_network.txt"), lang="julia").convert()
   ```

1. Add time-series data to [`experimental_data.jl`](https://github.com/pasmopy/breast_cancer/blob/master/training/erbb_network_jl/experimental_data.jl)

1. Set an objective function to be minimized in [`fitness.jl`](https://github.com/pasmopy/breast_cancer/blob/master/training/erbb_network_jl/fitness.jl)

1. Run [`optimize_parallel.sh`](https://github.com/pasmopy/breast_cancer/blob/master/training/optimize_parallel.sh)

   ```bash
   $ mv erbb_network_jl training
   $ cd training
   $ mkdir errout
   $ sh optimize_parallel.sh  # It will take more than a few days to optimize parameters.
   $ cd ..
   ```

   When finished, run:

   ```bash
   $ julia
   ```

   ```julia
   using BioMASS

   param2biomass("training")
   ```

   And you will get [`dat2npy/out/`](https://github.com/pasmopy/breast_cancer/tree/master/training/erbb_network_jl/dat2npy/out).
   This is the estimated parameter sets that [`biomass`](https://github.com/biomass-dev/biomass) can recognize and read.
   Copy `out/` to each patient-specific model folder via:

   ```python
   import os
   import shutil


   breast_cancer_models = []
   path_to_models = os.path.join("models", "breast")
   for model in os.listdir(path_to_models):
       if os.path.isdir(os.path.join(path_to_models, model)) and (
           model.startswith("TCGA_") or model.endswith("_BREAST")
       ):
           breast_cancer_models.append(model)
   # Set optimized parameters
   for model in breast_cancer_models:
       shutil.copytree(
           os.path.join("training", "erbb_network_jl", "dat2npy", "out"),
           os.path.join(path_to_models, f"{model}", "out"),
       )
   ```

## Patient-specific simulations

- Use `pasmopy.PatientModelSimulations`

  ```python
  import os
  import shutil
  from pathlib import Path

  from pasmopy import PatientModelSimulations

  import models.breast


  TCGA_ID = [
      l.strip() for l in Path("models", "breast", "sample_names.txt").read_text("utf-8").splitlines()
  ]
  # Create patient-specific models
  for patient in TCGA_ID:
      if patient != "TCGA_3C_AALK_01A":
          shutil.copytree(
              os.path.join("models", "breast", "TCGA_3C_AALK_01A"),
              os.path.join("models", "breast", f"{patient}"),
          )
  # Execute patient-specific models
  simulations = PatientModelSimulations(models.breast.__package__, TCGA_ID)
  simulations.run()
  ```