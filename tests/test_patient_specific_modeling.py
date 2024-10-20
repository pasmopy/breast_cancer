import os
import random
import shutil
import sys
import time
from typing import Callable, List, Optional

import numpy as np
from pasmopy import Model, PatientModelAnalyses, PatientModelSimulations, Text2Model
from pasmopy.preprocessing import WeightingFactors

from .C import *

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


with open(path_to_patient("selected_tnbc.txt"), mode="r") as f:
    TNBC_ID = f.read().splitlines()


def test_model_construction():
    # building a model
    try:
        shutil.rmtree(os.path.join("models", "erbb_network"))
    except FileNotFoundError:
        pass
    Text2Model(os.path.join("models", "erbb_network.txt")).convert()
    try:
        from models import erbb_network
    except ImportError:
        print("can't import erbb_network from models.")

    model = Model(erbb_network.__package__).create()
    # add weighting factors
    gene_expression = {
        "ErbB1": ["EGFR"],
        "ErbB2": ["ERBB2"],
        "ErbB3": ["ERBB3"],
        "ErbB4": ["ERBB4"],
        "Grb2": ["GRB2"],
        "Shc": ["SHC1", "SHC2", "SHC3", "SHC4"],
        "RasGAP": ["RASA1", "RASA2", "RASA3"],
        "PI3K": ["PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG"],
        "PTEN": ["PTEN"],
        "SOS": ["SOS1", "SOS2"],
        "Gab1": ["GAB1"],
        "RasGDP": ["HRAS", "KRAS", "NRAS"],
        "Raf": ["ARAF", "BRAF", "RAF1"],
        "MEK": ["MAP2K1", "MAP2K2"],
        "ERK": ["MAPK1", "MAPK3"],
        "Akt": ["AKT1", "AKT2"],
        "PTP1B": ["PTPN1"],
        "GSK3b": ["GSK3B"],
        "DUSP": ["DUSP5", "DUSP6", "DUSP7"],
        "cMyc": ["MYC"],
    }

    weighting_factors = WeightingFactors(model, gene_expression)
    weighting_factors.add_to_params()
    weighting_factors.set_search_bounds()
    # Edit observable.py, update normalization
    with open(
        os.path.join("models", "erbb_network", "observable.py"),
        mode="r",
        encoding="utf-8",
    ) as f:
        lines = f.readlines()
    for line_num, line in enumerate(lines):
        if line.startswith(f"{' ' * 8}self.normalization: dict = {{}}"):
            lines[line_num] = (
                f"{' ' * 8}self.normalization: dict = {{}}\n"
                + f"{' ' * 8}for obs_name in self.obs_names:\n"
                + f"{' ' * 12}self.normalization[obs_name] = "
                + f"{{'timepoint': None, 'condition': ['EGF', 'HRG']}}\n"
            )
    with open(
        os.path.join("models", "erbb_network", "observable.py"),
        mode="w",
        encoding="utf-8",
    ) as f:
        f.writelines(lines)
    # Edit set_search_param.py
    with open(
        os.path.join("models", "erbb_network", "set_search_param.py"),
        mode="r",
        encoding="utf-8",
    ) as f:
        lines = f.readlines()
    for line_num, line in enumerate(lines):
        if line.startswith("from .name2idx import C, V"):
            lines[line_num - 1] = f"{REQUIREMENTS}"
        elif line.startswith("class SearchParam(object):"):
            lines[line_num - 1] = f"\n{INDIVIDUALIZATION}\n\n"
        elif line.startswith(f"{' ' * 8}self.idx_initials = []"):
            lines[line_num] = f"{' ' * 8}self.idx_initials = [V.PIP2]\n"
        elif line.startswith(f"{' ' * 8}# parameter constraints"):
            lines[line_num - 1] = f"\n{INCORPORATION}\n"
    with open(
        os.path.join("models", "erbb_network", "set_search_param.py"),
        mode="w",
        encoding="utf-8",
    ) as f:
        f.writelines(lines)
    try:
        shutil.rmtree(os.path.join(PATH_TO_MODELS, "TCGA_3C_AALK_01A"))
    except FileNotFoundError:
        pass
    shutil.move(os.path.join("models", "erbb_network"), os.path.join(PATH_TO_MODELS, "TCGA_3C_AALK_01A"))
    assert os.path.isdir(path_to_patient("TCGA_3C_AALK_01A"))
    try:
        from models.breast import TCGA_3C_AALK_01A
    finally:
        model = Model(TCGA_3C_AALK_01A.__package__).create()
        # 220 parameters to be estimated & initial amount of PIP2.
        assert len(model.problem.idx_params) + len(model.problem.idx_initials) == 221


def test_patient_model_simulations(
    n_patients: int = 3,
    dynamical_feature: Optional[List[str]] = None,
):
    if not 1 <= n_patients <= 6:
        raise ValueError("`n_patients` must be lie within [1, 6].")
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
    simulations = PatientModelSimulations(
        models.breast.__package__,
        random.sample(TNBC_ID, n_patients),
    )
    start = time.time()
    assert simulations.run(n_proc=2) is None
    elapsed = time.time() - start
    print(f"Computation time for simulating {n_patients} patients: {elapsed/60:.1f} [min]")
    # Add new response characteristics
    get_droprate: Callable[[np.ndarray], float] = lambda time_course: -(time_course[-1] - np.max(time_course)) / (
        len(time_course) - np.argmax(time_course)
    )
    simulations.response_characteristics["droprate"] = get_droprate
    # Extract response characteristics and visualize patient classification
    if dynamical_feature is None:
        dynamical_feature = ["AUC", "droprate"]
    simulations.subtyping(
        "subtype_classification.pdf",
        {
            "Phosphorylated_Akt": {"EGF": dynamical_feature, "HRG": dynamical_feature},
            "Phosphorylated_ERK": {"EGF": dynamical_feature, "HRG": dynamical_feature},
            "Phosphorylated_c-Myc": {"EGF": dynamical_feature, "HRG": dynamical_feature},
        },
    )
    obs_names = ["Phosphorylated_Akt", "Phosphorylated_ERK", "Phosphorylated_c-Myc"]
    for observable in obs_names:
        assert os.path.isfile(os.path.join("classification", f"{observable}.csv"))
    assert os.path.isfile("subtype_classification.pdf")


def test_patient_model_analyses():
    for patient in TNBC_ID:
        assert os.path.isdir(os.path.join(PATH_TO_MODELS, patient))
    patients = random.sample(TNBC_ID, 1)
    analyses = PatientModelAnalyses(
        models.breast.__package__,
        patients,
        biomass_kws={"metric": "maximum", "style": "heatmap", "options": {"excluded_initials": ["PIP2"]}},
    )
    assert analyses.run(n_proc=2) is None
    for patient in patients:
        assert os.path.isfile(
            os.path.join(
                PATH_TO_MODELS,
                patient,
                "sensitivity_coefficients",
                "initial_condition",
                "maximum.npy",
            )
        )


def test_four_breast_cancer_cell_line_models():
    cell_lines: List[str] = []
    for model in os.listdir(os.path.join("models", "breast")):
        if model.endswith("_BREAST"):
            cell_lines.append(model)
    simulations = PatientModelSimulations(models.breast.__package__, cell_lines)
    assert simulations.run(n_proc=2) is None
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
        if patient in os.listdir(PATH_TO_MODELS):
            shutil.rmtree(path_to_patient(f"{patient}"))
    # patient classification
    files = os.listdir("classification")
    for file in files:
        if file.endswith(".csv"):
            os.remove(os.path.join("classification", f"{file}"))
    if os.path.isfile("subtype_classification.pdf"):
        os.remove("subtype_classification.pdf")
