import numpy as np
from biomass.dynamics.solver import solve_ode

from .name2idx import C, V
from .set_model import DifferentialEquation


class Observable(DifferentialEquation):
    """
    Correlating model simulations and experimental measurements.
    Attributes
    ----------
    obs_names : list of strings
        Names of model observables.
    t : range
        Simulation time span.
    conditions : list of strings
        Expetimental conditions.
    simulations : numpy.ndarray
        The numpy array to store simulation results.
    normalization : nested dict
        * 'timepoint' : Optional[int]
            The time point at which simulated values are normalized.
            If None, the maximum value will be used for normalization.
        * 'condition' : list of strings
            The experimental conditions to use for normalization.
            If empty, all conditions defined in sim.conditions will be used.
    experiments : list of dict
        Time series data.
    error_bars : list of dict
        Error bars to show in figures.
    """

    def __init__(self):
        super(Observable, self).__init__(perturbation={})
        self.obs_names: list = [
            "Phosphorylated_Akt",
            "Phosphorylated_ERK",
            "Phosphorylated_c-Myc",
        ]
        self.t: range = range(0, 120 + 1)
        self.conditions: list = [
            "EGF",
            "HRG",
        ]
        self.simulations: np.ndarray = np.empty(
            (len(self.obs_names), len(self.t), len(self.conditions))
        )
        self.normalization = {}
        for obs_name in self.obs_names:
            self.normalization[obs_name] = {"timepoint": None, "condition": ["EGF", "HRG"]}
        self.experiments: list = [None] * len(self.obs_names)
        self.error_bars: list = [None] * len(self.obs_names)

    def simulate(self, x, y0, _perturbation={}):
        if _perturbation:
            self.perturbation = _perturbation
        # Simulate model until all species reach steady state (if needed).
        # Define the untreated condition here.
        # y0 = self._get_steady_state(self.diffeq, y0, tuple(x))
        # if not y0:
        #    return False
        for i, condition in enumerate(self.conditions):
            if condition == "EGF":
                y0[V.EGF] = 10.0
                y0[V.HRG] = 0.0
            elif condition == "HRG":
                y0[V.EGF] = 0.0
                y0[V.HRG] = 10.0

            sol = solve_ode(self.diffeq, y0, self.t, tuple(x))

            if sol is None:
                return False
            else:
                self.simulations[self.obs_names.index("Phosphorylated_Akt"), :, i] = sol.y[
                    V.AktP, :
                ]
                self.simulations[self.obs_names.index("Phosphorylated_ERK"), :, i] = (
                    sol.y[V.ERKP, :] + sol.y[V.ERKPP, :]
                )
                self.simulations[self.obs_names.index("Phosphorylated_c-Myc"), :, i] = (
                    sol.y[V.cMycP, :] + sol.y[V.cMycPP, :]
                )

    def set_data(self):
        self.experiments[self.obs_names.index("Phosphorylated_Akt")] = {
            "EGF": [0.0, 1.0, 0.409, 0.229, 0.057, 0.03, 0.005, 0.0],
            "HRG": [0.0, 0.017, 0.084, 0.03, 0.023, 0.038, 0.028, 0.024],
        }
        self.error_bars[self.obs_names.index("Phosphorylated_Akt")] = {
            "EGF": [0.017, 0.0, 0.133, 0.077, 0.026, 0.018, 0.015, 0.032],
            "HRG": [0.017, 0.015, 0.05, 0.018, 0.014, 0.017, 0.008, 0.031],
        }

        self.experiments[self.obs_names.index("Phosphorylated_ERK")] = {
            "EGF": [0.278, 0.884, 1.0, 0.993, 0.817, 0.672, 0.419, 0.369],
            "HRG": [0.278, 0.077, 0.182, 0.157, 0.0, 0.007, 0.066, 0.293],
        }
        self.error_bars[self.obs_names.index("Phosphorylated_ERK")] = {
            "EGF": [0.051, 0.095, 0.054, 0.081, 0.07, 0.088, 0.064, 0.027],
            "HRG": [0.051, 0.045, 0.067, 0.06, 0.024, 0.049, 0.079, 0.149],
        }

        self.experiments[self.obs_names.index("Phosphorylated_c-Myc")] = {
            "EGF": [0.229, 0.235, 0.0, 0.461, 0.814, 0.969, 1.0, 0.928],
            "HRG": [0.229, 0.167, 0.187, 0.233, 0.361, 0.539, 0.719, 0.841],
        }
        self.error_bars[self.obs_names.index("Phosphorylated_c-Myc")] = {
            "EGF": [0.285, 0.361, 0.349, 0.38, 0.206, 0.188, 0.186, 0.192],
            "HRG": [0.285, 0.091, 0.095, 0.05, 0.078, 0.145, 0.206, 0.077],
        }

    def get_timepoint(self, obs_name):
        if obs_name in self.obs_names:
            return [0, 5, 15, 30, 45, 60, 90, 120]
