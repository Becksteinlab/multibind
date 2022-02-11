from multibind.multibind import Multibind
import pandas as pd
import tempfile
import os
import numpy as np
from math import isclose


def rate_matrix(filename: str) -> (Multibind, np.ndarray, np.ndarray):
    """Build a rate matrix from a connectivity file with rates.

    The rates are converted into free energy differences between states and
    then made thermodynamically consistent.

    Parameters
    ----------
    filename : str
        path to the CSV file containing the state names and connections in the form of rates
        Note that the order of input rates matters.

    Returns
    -------
    graph : Multibind
        Multibind graph with the appropriate states and connections filled in from the input rate file
    rate_matrix : 2D array
        2D array showing the rate of going from the state in index 0 to the state in index 1.
    rate_std_err_matrix : 2D array
        2D array with the standard error on the projected rates.
    """

    _rates = pd.read_csv(filename, header=0)

    with tempfile.TemporaryDirectory() as tmpdirname:
        # multibind requires that states and graphs are defined within csv files, we will create temporary files
        states_filename = os.path.join(tmpdirname, "states.csv")
        graph_filename = os.path.join(tmpdirname, "graph.csv")

        _states = []
        _edges = []

        for index, row in _rates.iterrows():
            state1 = row['state1']
            state2 = row['state2']

            # find all unique states
            _states.append(state1) if state1 not in _states else None
            _states.append(state2) if state2 not in _states else None

            # find all unique edges. Implicitly assumes that forward rates appear first...
            # todo add this to the docstring
            _edges.append((state1, state2)) if ((state1, state2) not in _edges and (state2, state1) not in _edges) else None

        states = pd.DataFrame({'name': _states})  # create dataframe for states, which only needs the names

        # lists of values that will be added to the graph dataframe
        state1_list = []
        state2_list = []
        value_list = []
        variance_list = []
        ligand_list = []
        standard_state_list = []

        for state1, state2 in _edges:
            state1_list.append(state1)
            state2_list.append(state2)

            k_f, var_f = _rates[(_rates.state1 == state1) & (_rates.state2 == state2)].values[0, 2:]
            k_b, var_b = _rates[(_rates.state1 == state2) & (_rates.state2 == state1)].values[0, 2:]

            dg = np.log(k_b / k_f)
            var = var_f / k_f ** 2 + var_b / k_b ** 2

            value_list.append(dg)
            variance_list.append(var)
            ligand_list.append("helm")
            standard_state_list.append(1)

            graph = pd.DataFrame({"state1": state1_list,
                                  "state2": state2_list,
                                  "value": value_list,
                                  "variance": variance_list,
                                  "ligand": ligand_list,
                                  "standard_state": standard_state_list})

        states.to_csv(states_filename, index=False)
        graph.to_csv(graph_filename, index=False)

        c = Multibind()
        c.read_graph(graph_filename)
        c.read_states(states_filename)

        c.build_cycle(pH=0)
        c.MLE()

        states, connections, g = c.states.name.values, c.graph[['state1', 'state2']].values, c.g_mle
        standard_errors = np.sqrt(np.diagonal(c.covariance_matrix))

        n = states.shape[0]
        _rate_matrix = np.zeros((n, n))
        _rate_matrix_SE = np.zeros_like(_rate_matrix)

        for i, j in connections:
            # indices of states
            _i = np.argwhere(states == i)[0, 0]
            _j = np.argwhere(states == j)[0, 0]

            # input rates and their variances
            kijb, vij = _rates[(_rates.state1 == i) & (_rates.state2 == j)].values[0, 2:]
            kjib, vji = _rates[(_rates.state1 == j) & (_rates.state2 == i)].values[0, 2:]

            # standard errors output by multibind
            si = standard_errors[_i]
            sj = standard_errors[_j]

            # most functions will take in standard error instead of variance
            sijb = np.sqrt(vij)
            sjib = np.sqrt(vji)

            sijb, sjib = clean_sigmas(sijb, sjib)

            # free energy difference
            dg = g[_j] - g[_i]

            # standard errors of forward and backward rates
            skji = kji_standard_error(kijb, kjib, sijb, sjib, dg, sj, si)
            skij = kij_standard_error(kijb, kjib, sijb, sjib, dg, sj, si)

            # new rates
            kji, kij = project_rates(kijb, kjib, sijb, sjib, dg)

            # confirm that we can recreate the target free energy differenc
            assert np.isclose(dg, np.log(kji / kij)), (dg, np.log(kji / kij))

            _rate_matrix[_i, _j] = kij
            _rate_matrix[_j, _i] = kji
            _rate_matrix_SE[_i, _j] = skij
            _rate_matrix_SE[_j, _i] = skji

        return c, _rate_matrix, _rate_matrix_SE


def clean_sigmas(r1: float, r2: float) -> (float, float):
    """If sigmas are zero, it's best to just make them close to zero to avoid errors.

    Parameters
    ----------
    r1 : float
        First rate
    r2 : float
        Second rate

    Returns
    -------
    (float, float)
        The potentially transformed rates
    """
    if isclose(r1, 0):
        r1 = 1e-12
    if isclose(r2, 0):
        r2 = 1e-12
    return r1, r2


def project_rates(kijb: float, kjib: float, sijb: float, sjib: float, dg: float, infinity: float = 1e12) -> (float, float):
    """Project input rates and standard errors onto the consistency line.

    Parameters
    ----------
    kijb : float
        Input rate from state i to j
    kjib : float
        Input rate from state j to i
    sijb : float
        Input standard error the rate from state i to j
    sjib : float
        Input standard error the rate from state j to i
    dg : float
        Target free energy that defines the consistency line
    infinity: float, optional
        Value to be used for s_ijji or s_jiij if any input variances are zero

    Returns
    -------
    (kji, kij) : (float, float)
        The new projected rates
    """
    b_weight = np.exp(-dg)

    vij = sijb**2
    vji = sjib**2

    s_ijji = vij / vji
    s_jiij = 1 / s_ijji

    kji = kjib / (1 + s_jiij * b_weight**2) + kijb / (s_ijji / b_weight + b_weight)
    kij = kji * b_weight

    return (kji, kij)


def dG_standard_error(sj: float, si: float) -> float:
    """The standard error of the free energy difference between states i and j.
    Note that this is a symmetric function.

    Parameters
    ----------
    sj : float
        Standard error for the estimate of the free energy of state j
    si : float
        Standard error for the estimate of the free energy of state i

    Returns
    -------
    sdg : float
        Standard error of the free energy difference
    """
    return np.sqrt(sj**2 + si**2)


def kji_standard_error(kijb: float, kjib: float, sijb: float, sjib: float, target_dG: float, sj: float, si: float) -> float:
    """The standard error of kji from the projection method.

    Parameters
    ----------
    kijb : float
        Input rate from state i to j
    kjib : float
        Input rate from state j to i
    sijb : float
        Input standard error for the rate from state i to j
    sjib : float
        Input standard error for the rate from state j to i
    target_dG : float
        The determined free energy difference that defines the consistency line.
    sj : float
        Standard error of the free energy of state j
    si : float
        Standard error of the free energy of state i

    Returns
    -------
    sji : float
        Standard error of the projected kji rate
    """
    vji = sjib**2
    vij = sijb**2
    b_weight = np.exp(-target_dG)
    sdg = np.sqrt(sj**2 + si**2)

    skji = np.sqrt((b_weight**-2 * vji**2 * sdg**2 * (-b_weight**-2 * kijb * vij + 2 * b_weight**-1 * kjib * vij + kijb * vji)**2) / (b_weight**-2 * vij + vji)**4)

    return skji


def kij_standard_error(kijb: float, kjib: float, sijb: float, sjib: float, target_dG: float, sj: float, si: float) -> float:
    """The standard error of kij from the projection method.

    Parameters
    ----------
    kijb : float
        Input rate from state i to j
    kjib : float
        Input rate from state j to i
    sijb : float
        Input standard error for the rate from state i to j
    sjib : float
        Input standard error for the rate from state j to i
    target_dG : float
        The determined free energy difference that defines the consistency line.
    sj : float
        Standard error of the free energy of state j
    si : float
        Standard error of the free energy of state i

    Returns
    -------
    sij : float
        Standard error of the projected kij rate
    """
    sdg = dG_standard_error(sj, si)
    b_weight = np.exp(-target_dG)

    skji = kji_standard_error(kijb, kjib, sijb, sjib, target_dG, sj, si)

    _k_ji, _k_ij = project_rates(kijb, kjib, sijb, sjib, target_dG)

    skij = np.sqrt(sdg**2 * (-_k_ji * b_weight)**2
                   + skji**2 * b_weight**2)
    return skij
