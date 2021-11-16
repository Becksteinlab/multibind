from multibind.multibind import Multibind
import pandas as pd
import tempfile
import os
import numpy as np


def rate_matrix(filename: str):
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

        n = states.shape[0]
        _rate_matrix = np.zeros((n, n))

        for i, j in connections:
            _i = np.argwhere(states == i)[0, 0]
            _j = np.argwhere(states == j)[0, 0]
            k_ij, var_ij = _rates[(_rates.state1 == i) & (_rates.state2 == j)].values[0, 2:]
            k_ji, var_ji = _rates[(_rates.state1 == j) & (_rates.state2 == i)].values[0, 2:]
            dg = g[_j] - g[_i]
            b_weight = np.exp(-dg)
            try:
                s_ijji = var_ij / var_ji
            except ZeroDivisionError:
                s_ijji = 1e12
            try:
                s_jiij = 1 / s_ijji
            except ZeroDivisionError:
                s_jiij = 1e12

            _k_ji = k_ij / (b_weight + s_ijji / b_weight ** 2) + k_ji / (b_weight ** 2 * s_jiij + 1)
            _k_ij = _k_ji * b_weight

            assert np.isclose(dg, np.log(_k_ji / _k_ij)), (dg, np.log(_k_ji / _k_ij))

            _rate_matrix[_i, _j] = _k_ij
            _rate_matrix[_j, _i] = _k_ji

        return c, _rate_matrix
