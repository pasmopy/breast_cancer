---
layout: default
title: Subtype classification based on the ErbB signaling dynamics
parent: Getting started
nav_order: 4
---

# Subtype classification based on the ErbB signaling dynamics

## Extraction of response characteristics from patient-specific simulations

- Execute `subtyping()`

  ```python
  simulations.subtyping(
      fname=None,
      dynamical_features={
          "Phosphorylated_Akt": {"EGF": ["max"], "HRG": ["max"]},
          "Phosphorylated_ERK": {"EGF": ["max"], "HRG": ["max"]},
          "Phosphorylated_c-Myc": {"EGF": ["max"], "HRG": ["max"]},
      }
  )
  ```

## Model-based patient stratification

- Run [`brca_heatmap.R`](https://github.com/pasmopy/breast_cancer/tree/master/classification/brca_heatmap.R)

  ```bash
  $ cd classification
  # $ Rscript brca_heatmap.R [n_cluster: int] [figsize: tuple]
  $ Rscript brca_heatmap.R 6 8,5
  ```