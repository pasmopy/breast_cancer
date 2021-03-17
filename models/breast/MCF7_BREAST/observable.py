from typing import Callable, List, Optional

import numpy as np
from scipy.integrate import solve_ivp

from .name2idx import C, V
from .set_model import DifferentialEquation

observables = [
    "Phosphorylated_Akt",
    "Phosphorylated_ERK",
    "Phosphorylated_cMyc",
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
                self.simulations[observables.index("Phosphorylated_ERK"), :, i] = sol.y[V.ERKP, :] + sol.y[V.ERKPP, :]
                self.simulations[observables.index("Phosphorylated_cMyc"), :, i] = (
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
        data_pAkt = self._norm01(
            egf1=[32101, 156970, 90301, 76709, 63640, 52536, 46414, 57329],
            hrg1=[32101, 565508, 551901, 560064, 489678, 408802, 425323, 451502],
            egf2=[11612, 96189, 43622, 43238, 41007, 29902, 19255, 35079],
            hrg2=[11612, 397931, 432609, 417622, 434519, 509919, 361041, 292523],
            egf3=[66038, 208525, 102689, 117308, 125158, 92086, 68587, 78252],
            hrg3=[66038, 563079, 573540, 521062, 447462, 383774, 434807, 409615],
        )
        self.experiments[observables.index("Phosphorylated_Akt")] = {
            "EGF": data_pAkt[0],
            "HRG": data_pAkt[1],
        }
        self.error_bars[observables.index("Phosphorylated_Akt")] = {
            "EGF": data_pAkt[2],
            "HRG": data_pAkt[3],
        }

        data_pERK = self._norm01(
            egf1=[65481, 446949, 221435, 283171, 265152, 266056, 204912, 188972],
            hrg1=[65481, 698717, 766252, 710005, 693622, 691856, 522173, 334410],
            egf2=[41927, 507623, 169918, 193671, 164088, 145916, 110844, 130362],
            hrg2=[41927, 605118, 699511, 654697, 579863, 490649, 299946, 229297],
            egf3=[118995, 807929, 338665, 267160, 253820, 230200, 157620, 153112],
            hrg3=[118995, 710436, 673318, 615206, 612686, 523198, 390301, 257664],
        )
        self.experiments[observables.index("Phosphorylated_ERK")] = {
            "EGF": data_pERK[0],
            "HRG": data_pERK[1],
        }
        self.error_bars[observables.index("Phosphorylated_ERK")] = {
            "EGF": data_pERK[2],
            "HRG": data_pERK[3],
        }
        '''
        data_pcFos = self._norm01(
            egf1=[43200, 101848, 134050, 187681, 274701, 188891, 186912, 147868],
            hrg1=[43200, 243299, 340259, 537344, 583257, 527613, 551327, 630883],
            egf2=[36344, 99849, 173325, 179897, 207943, 155466, 154118, 138196],
            hrg2=[36344, 139813, 245333, 389460, 402734, 556006, 513591, 432916],
            egf3=[43604, 83374, 108733, 116103, 113879, 95504, 94969, 94662],
            hrg3=[43604, 111136, 343365, 464180, 453578, 440094, 404483, 589354],
        )
        self.experiments[observables.index("Phosphorylated_cFos")] = {
            "EGF": data_pcFos[0],
            "HRG": data_pcFos[1],
        }
        self.error_bars[observables.index("Phosphorylated_cFos")] = {
            "EGF": data_pcFos[2],
            "HRG": data_pcFos[3],
        }
        '''
        data_pcMyc = self._norm01(
            egf1=[115975, 226001, 166894, 194150, 263331, 235172, 126949, 91142],
            hrg1=[115975, 62515, 81364, 155844, 390689, 664641, 848356, 856941],
            egf2=[185069, 276202, 204012, 234391, 290020, 360762, 325531, 242455],
            hrg2=[185069, 234416, 251732, 333993, 550670, 859790, 939956, 769616],
            egf3=[127244, 186118, 163387, 132053, 192949, 220987, 184381, 151547],
            hrg3=[127244, 110676, 152880, 277206, 461217, 637033, 908235, 712427],
        )
        self.experiments[observables.index("Phosphorylated_cMyc")] = {
            "EGF": data_pcMyc[0],
            "HRG": data_pcMyc[1],
        }
        self.error_bars[observables.index("Phosphorylated_cMyc")] = {
            "EGF": data_pcMyc[2],
            "HRG": data_pcMyc[3],
        }

    @staticmethod
    def get_timepoint(obs_name):
        if obs_name in observables:
            return [0, 5, 15, 30, 45, 60, 90, 120]
