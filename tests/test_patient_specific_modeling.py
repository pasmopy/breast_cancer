import os
import random
import shutil
import sys
import time
from typing import Callable, List

import numpy as np
from pasmopy import PatientModelSimulations

try:
    import models.breast
except ImportError:
    print("can't import 'breast' from 'models'.")

if sys.version_info[:2] < (3, 7):
    raise RuntimeError("`pasmopy` requires Python 3.7+ to run.")

PATH_TO_MODELS: str = os.path.join("models", "breast")


def path_to_patient(patient_id: str) -> str:
    return os.path.join(PATH_TO_MODELS, patient_id)


with open(path_to_patient("sample_names.txt"), mode="r") as f:
    TCGA_ID = f.read().splitlines()


def test_patient_model_simulations():
    # Initialization
    for patient in TCGA_ID:
        if patient in os.listdir(PATH_TO_MODELS) and patient != "TCGA_3C_AALK_01A":
            shutil.rmtree(path_to_patient(f"{patient}"))
    try:
        shutil.rmtree(os.path.join(path_to_patient("TCGA_3C_AALK_01A"), "out"))
    except FileNotFoundError:
        pass
    # Set optimized parameter sets
    breast_cancer_models: List[str] = []
    for f in os.listdir(PATH_TO_MODELS):
        if os.path.isdir(path_to_patient(f)) and (f.startswith("TCGA_") or f.endswith("_BREAST")):
            breast_cancer_models.append(f)
    for model in breast_cancer_models:
        if os.path.isdir(os.path.join(path_to_patient(f"{model}"), "out")):
            shutil.rmtree(os.path.join(path_to_patient(f"{model}"), "out"))
        shutil.copytree(
            os.path.join("training", "erbb_network_jl", "dat2npy", "out"),
            os.path.join(path_to_patient(f"{model}"), "out"),
        )
    # Create patient-specific models
    for patient in TCGA_ID:
        if patient != "TCGA_3C_AALK_01A":
            shutil.copytree(path_to_patient("TCGA_3C_AALK_01A"), path_to_patient(f"{patient}"))
    # Execute patient-specific models
    simulations = PatientModelSimulations(models.breast.__package__, random.sample(TCGA_ID, 10))
    start = time.time()
    assert simulations.run() is None
    elapsed = time.time() - start
    print(f"Computation time for simulating 10 patients: {elapsed/60:.1f} [min]")
    # Add new response characteristics
    _droprate: Callable[[np.ndarray], float] = lambda time_course: -(
        time_course[-1] - np.max(time_course)
    ) / (len(time_course) - np.argmax(time_course))
    simulations.response_characteristics["droprate"] = _droprate
    simulations.response_characteristics["argmax"] = np.argmax
    # Extract response characteristics and visualize patient classification
    simulations.subtyping(
        "subtype_classification.pdf",
        {
            "Phosphorylated_Akt": {"EGF": ["max"], "HRG": ["AUC"]},
            "Phosphorylated_ERK": {"EGF": ["droprate"], "HRG": ["droprate"]},
            "Phosphorylated_c-Myc": {"EGF": ["argmax"], "HRG": ["argmax"]},
        },
    )
    observables = ["Phosphorylated_Akt", "Phosphorylated_ERK", "Phosphorylated_c-Myc"]
    for obs in observables:
        assert os.path.isfile(os.path.join("classification", f"{obs}.csv"))
    assert os.path.isfile("subtype_classification.pdf")


def test_four_breast_cancer_cell_line_models():
    cell_lines: List[str] = []
    for model in os.listdir(os.path.join("models", "breast")):
        if model.endswith("_BREAST"):
            cell_lines.append(model)
    simulations = PatientModelSimulations(models.breast.__package__, cell_lines)
    assert simulations.run() is None
    for model in cell_lines:
        simulated_values = np.load(
            os.path.join(
                "models",
                "breast",
                model,
                "simulation_data",
                "simulations_all.npy",
            )
        )
        assert np.isfinite(simulated_values).all()


def test_cleanup_models():
    # patients
    for patient in TCGA_ID:
        if patient in os.listdir(PATH_TO_MODELS) and patient != "TCGA_3C_AALK_01A":
            shutil.rmtree(path_to_patient(f"{patient}"))
    # parameter sets
    shutil.rmtree(os.path.join(path_to_patient("TCGA_3C_AALK_01A"), "out"))
    # patient classification
    files = os.listdir("classification")
    for file in files:
        if file.endswith(".csv"):
            os.remove(os.path.join("classification", f"{file}"))
    os.remove("subtype_classification.pdf")
