from typing import Callable, List, Optional

import numpy as np
from scipy.integrate import solve_ivp

from .name2idx import C, V
from .set_model import DifferentialEquation

observables = [
    "Phosphorylated_Akt",
    "Phosphorylated_ERK",
    "Phosphorylated_c-Myc",
]


class NumericalSimulation(DifferentialEquation):
    """Simulate a model using scipy.integrate.ode

    Attributes
    ----------
    normalization : nested dict
        Keys for each observable
        ------------------------
        * 'timepoint' : Optional[int]
            The time point at which simulated values are normalized.
            If None, the maximum value will be used for normalization.

        * 'condition' : list of strings
            The experimental conditions to use for normalization.
            If empty, all conditions defined in sim.conditions will be used.

    """

    def __init__(self):
        super().__init__(perturbation={})
        self.normalization = {}
        for obs_name in observables:
            self.normalization[obs_name] = {"timepoint": None, "condition": []}

    t = range(0, 120 + 1)
    # Experimental conditions
    conditions = [
        "EGF",
        "HRG",
    ]

    simulations = np.empty((len(observables), len(t), len(conditions)))

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

            sol = self._solveode(self.diffeq, y0, self.t, tuple(x))

            if sol is None:
                return False
            else:
                self.simulations[observables.index("Phosphorylated_Akt"), :, i] = sol.y[V.AktP, :]
                self.simulations[observables.index("Phosphorylated_ERK"), :, i] = (
                    sol.y[V.ERKP, :] + sol.y[V.ERKPP, :]
                )
                self.simulations[observables.index("Phosphorylated_c-Myc"), :, i] = (
                    sol.y[V.cMycP, :] + sol.y[V.cMycPP, :]
                )

    @staticmethod
    def _solveode(
        diffeq: Callable,
        y0: List[float],
        t: range,
        f_params: tuple,
        method: str = "BDF",
        options: Optional[dict] = None,
    ):
        """
        Solve a system of ordinary differential equations.
        Parameters
        ----------
        diffeq : callable f(t, y, *x)
            Right-hand side of the differential equation.
        y0 : array
            Initial condition on y (can be a vector).
        t : array
            A sequence of time points for which to solve for y.
        f_params : tuple
            Model parameters.
        method : str (default: "BDF")
            Integration method to use.
        options : dict, optional
            Options passed to a chosen solver.
        Returns
        -------
        sol : OdeResult
            Represents the solution of ODE.
        """
        if options is None:
            options = {}
        options.setdefault("rtol", 1e-8)
        options.setdefault("atol", 1e-8)
        try:
            sol = solve_ivp(
                diffeq,
                (t[0], t[-1]),
                y0,
                method=method,
                t_eval=t,
                args=f_params,
                **options,
            )
            return sol if sol.success else None
        except ValueError:
            return None

    def _get_steady_state(
        self,
        diffeq: Callable,
        y0: List[float],
        f_params: tuple,
        eps: float = 1e-6,
    ):
        """
        Find the steady state for the untreated condition.
        Parameters
        ----------
        diffeq : callable f(t, y, *x)
            Right-hand side of the differential equation.
        y0 : array
            Initial condition on y (can be a vector).
        f_params : tuple
            Model parameters.
        eps : float (default: 1e-6)
            Run until a time t for which the maximal absolute value of the
            regularized relative derivative was smaller than eps.
        Returns
        -------
        y0 : array
            Steady state concentrations of all species.
        """
        while True:
            sol = self._solveode(diffeq, y0, range(2), f_params)
            if sol is None or np.max(np.abs((sol.y[:, -1] - y0) / (np.array(y0) + eps))) < eps:
                break
            else:
                y0 = sol.y[:, -1].tolist()

        return [] if sol is None else sol.y[:, -1].tolist()


class ExperimentalData(object):
    """
    Set experimental data.

    Attributes
    ----------
    experiments : list of dict
        Time series data.

    error_bars : list of dict
        Error bars to show in figures.

    """

    def __init__(self):
        self.experiments = [None] * len(observables)
        self.error_bars = [None] * len(observables)

    @staticmethod
    def _norm01(egf1, hrg1, egf2, hrg2, egf3, hrg3):
        """Normalize data from 0 to 1 and Get standard error of the mean

        Parameters
        ----------
        egf%d : list
            EGF stimulation #%d
        hrg%d : list
            HRG stimulation #%d

        Returns
        -------
        egf_ave, hrg_ave, egf_sem, hrg_sem : tuple
            Averaged values and their standard error of the mean
        """
        data1 = np.array([egf1, hrg1])
        data2 = np.array([egf2, hrg2])
        data3 = np.array([egf3, hrg3])

        # max -> 1 to compare biological replicates
        data1 = data1 / np.max(data1)
        data2 = data2 / np.max(data2)
        data3 = data3 / np.max(data3)

        egf_ave = np.mean(np.stack([data1[0], data2[0], data3[0]]), axis=0)
        hrg_ave = np.mean(np.stack([data1[1], data2[1], data3[1]]), axis=0)

        ave_vec = np.stack([egf_ave, hrg_ave])
        ave_min = np.min(ave_vec)
        ave_max = np.max(ave_vec)

        # To normalize max -> 1, min -> 0
        data1 = (data1 - ave_min) / (ave_max - ave_min)
        data2 = (data2 - ave_min) / (ave_max - ave_min)
        data3 = (data3 - ave_min) / (ave_max - ave_min)

        egf_ave = np.mean(np.stack([data1[0], data2[0], data3[0]]), axis=0)
        hrg_ave = np.mean(np.stack([data1[1], data2[1], data3[1]]), axis=0)
        egf_sem = np.std(np.stack([data1[0], data2[0], data3[0]]), axis=0, ddof=1) / (3 ** 0.5)
        hrg_sem = np.std(np.stack([data1[1], data2[1], data3[1]]), axis=0, ddof=1) / (3 ** 0.5)

        return egf_ave, hrg_ave, egf_sem, hrg_sem

    def set_data(self):
        pass

    @staticmethod
    def get_timepoint(obs_name):
        if obs_name in observables:
            return [0, 5, 15, 30, 45, 60, 90, 120]
