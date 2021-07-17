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
            "EGF": [0.0, 0.242, 0.087, 0.088, 0.082, 0.045, 0.017, 0.043],
            "HRG": [0.0, 0.976, 1.0, 0.96, 0.876, 0.836, 0.77, 0.719],
        }
        self.error_bars[self.obs_names.index("Phosphorylated_Akt")] = {
            "EGF": [0.031, 0.058, 0.033, 0.04, 0.047, 0.034, 0.027, 0.022],
            "HRG": [0.031, 0.08, 0.054, 0.056, 0.03, 0.117, 0.018, 0.075],
        }

        self.experiments[self.obs_names.index("Phosphorylated_ERK")] = {
            "EGF": [0.0, 0.794, 0.259, 0.269, 0.237, 0.216, 0.129, 0.13],
            "HRG": [0.0, 0.93, 1.0, 0.918, 0.866, 0.771, 0.512, 0.311],
        }
        self.error_bars[self.obs_names.index("Phosphorylated_ERK")] = {
            "EGF": [0.031, 0.144, 0.062, 0.032, 0.039, 0.047, 0.038, 0.023],
            "HRG": [0.031, 0.016, 0.066, 0.067, 0.05, 0.092, 0.091, 0.045],
        }

        self.experiments[self.obs_names.index("Phosphorylated_c-Myc")] = {
            "EGF": [0.011, 0.125, 0.058, 0.07, 0.151, 0.18, 0.099, 0.034],
            "HRG": [0.011, 0.0, 0.034, 0.156, 0.434, 0.765, 1.0, 0.848],
        }
        self.error_bars[self.obs_names.index("Phosphorylated_c-Myc")] = {
            "EGF": [0.023, 0.031, 0.013, 0.037, 0.038, 0.05, 0.07, 0.052],
            "HRG": [0.023, 0.062, 0.059, 0.061, 0.044, 0.074, 0.004, 0.079],
        }

    def get_timepoint(self, obs_name):
        if obs_name in self.obs_names:
            return [0, 5, 15, 30, 45, 60, 90, 120]
