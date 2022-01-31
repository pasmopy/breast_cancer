import os
import random
import shutil
import sys
import time
from typing import List

from pasmopy import PatientModelSimulations

try:
    import models.colon
except ImportError:
    print("can't import 'colon' from 'models'.")

if sys.version_info[:2] < (3, 7):
    raise RuntimeError("`pasmopy` requires Python 3.7+ to run.")

PATH_TO_MODELS: str = os.path.join("models", "colon")


def path_to_patient(patient_id: str) -> str:
    return os.path.join(PATH_TO_MODELS, patient_id)


with open(path_to_patient("sample_names.txt"), mode="r") as f:
    TCGA_ID = f.read().splitlines()


def test_patient_model_simulations():
    # Initialization
    for patient in TCGA_ID:
        if patient in os.listdir(PATH_TO_MODELS) and patient != "TCGA_4T_AA8H_01A":
            shutil.rmtree(path_to_patient(f"{patient}"))
    try:
        shutil.rmtree(os.path.join(path_to_patient("TCGA_4T_AA8H_01A"), "out"))
    except FileNotFoundError:
        pass
    # Set optimized parameter sets
    colon_cancer_models: List[str] = []
    for f in os.listdir(PATH_TO_MODELS):
        if os.path.isdir(path_to_patient(f)) and f.startswith("TCGA_"):
            colon_cancer_models.append(f)
    for model in colon_cancer_models:
        if os.path.isdir(os.path.join(path_to_patient(f"{model}"), "out")):
            shutil.rmtree(os.path.join(path_to_patient(f"{model}"), "out"))
        shutil.copytree(
            os.path.join("training", "erbb_network_jl", "dat2npy", "out"),
            os.path.join(path_to_patient(f"{model}"), "out"),
        )
    # Create patient-specific models
    for patient in TCGA_ID:
        if patient != "TCGA_4T_AA8H_01A":
            shutil.copytree(path_to_patient("TCGA_4T_AA8H_01A"), path_to_patient(f"{patient}"))
    # Execute patient-specific models
    simulations = PatientModelSimulations(models.colon.__package__, random.sample(TCGA_ID, 3))
    start = time.time()
    assert simulations.run() is None
    elapsed = time.time() - start
    print(f"Computation time for simulating 3 patients: {elapsed/60:.1f} [min]")


def test_cleanup_models():
    # patients
    for patient in TCGA_ID:
        if patient in os.listdir(PATH_TO_MODELS) and patient != "TCGA_4T_AA8H_01A":
            shutil.rmtree(path_to_patient(f"{patient}"))
    # parameter sets
    shutil.rmtree(os.path.join(path_to_patient("TCGA_4T_AA8H_01A"), "out"))
