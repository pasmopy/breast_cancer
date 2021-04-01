import os
import random
import shutil
import sys
import time
from typing import List

from dyaus_dev import PatientModelSimulations

if sys.version_info[:2] < (3, 7):
    raise RuntimeError("Python version >= 3.7 required.")

PATH_TO_MODELS: str = os.path.join("models", "breast")


def path_to_patient(path: str) -> str:
    return os.path.join(PATH_TO_MODELS, path)


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
    models: List[str] = []
    for f in os.listdir(PATH_TO_MODELS):
        if os.path.isdir(path_to_patient(f)) and (f.startswith("TCGA_") or f.endswith("_BREAST")):
            models.append(f)
    for model in models:
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
    simulations = PatientModelSimulations("models.breast", random.sample(TCGA_ID, 10))
    start = time.time()
    assert simulations.run() is None
    elapsed = time.time() - start
    print(f"Computation time for simulating 10 patients: {elapsed/60:.1f} [min]")


def test_cleanup_models():
    # patients
    for patient in TCGA_ID:
        if patient in os.listdir(PATH_TO_MODELS) and patient != "TCGA_3C_AALK_01A":
            shutil.rmtree(path_to_patient(f"{patient}"))
    # parameter sets
    shutil.rmtree(os.path.join(path_to_patient("TCGA_3C_AALK_01A"), "out"))
