import pandas as pd
import numpy as np
from numpy.linalg import pinv
import networkx as nx
from xarray import Dataset
from itertools import product
from collections import OrderedDict
from multibind.chem import protonation_free_energy, protonation_free_energy_standard_error, binding_free_energy_general


class InvalidConcentrationError(Exception):
    """Indicates that an invalid concentration was found.
    """
    pass


class MultibindScanner(object):
    """Run multibind across a range of concentrations.

    Attributes
    ----------
    results.deltas : xr.DataArray
        DataArray containing the input free energy differences across multiple ligand concentrations.
        See docstring for `Multibind` for more info.
    results.covariance : xr.DataArray
        DataArray containing the covariance between calculated state free energies across multiple ligand
        concentrations. See docstring for `Multibind` for more info.
    results.std_errors : xr.DataArray
        DataArray containing the standard error for each state free energy across multiple ligand concentrations.
        See docstring for `Multibind` for more info.
    """

    def __init__(self, statefile: str, graphfile: str, comment_char: str = '#'):
        self.c = Multibind(states_filename=statefile, graph_filename=graphfile, comment_char=comment_char)

    def run(self, concentrations: dict, svd: bool = True) -> None:
        """Create an xarray dataset containing thermodynamic properties as a
        function of varying concentrations.

        Parameters
        ----------
        concentrations : dict
            Dictionary containing the concentrations to scan over.

            Two input forms are valid. A dictionary with two items will be interpreted as an outer product.
            A dictionary with one item will be interpreted as specific points (not implemented yet).

            e.g. (('H+', 'Na+'): [(1, 0.250), (2, 0.150), (3, 0.100)]) will indicate that three graphs will be
            solved. One with pH = 1 and [Na+] = 0.250 M, another with pH = 2 and [Na+] = 0.150 M, and the last being pH = 3 and [Na+] = 0.100 M.

            {'H+' : [1, 2, 3, 4, 5], 'Na+' : [0.250, 0.150, 0.100]} will indicate that 15 graphs will be solved.
            The pairs will be the outer product of the two arrays.

        svd : bool (optional)
            Whether or not to use SVD to solve the graph. Otherwise use Newton-Raphson. Defaults to True.

        Returns
        -------
        xarray DataSet
        """

        if type(concentrations) is list or type(concentrations) is tuple:
            raise NotImplementedError

        elif isinstance(concentrations, dict):
            names = []
            values = []

            for k, v in concentrations.items():
                if k == 'H+':
                    names.append('pH')
                else:
                    names.append(k)
                values.append(v)

            grid = list(product(*values))
            indices = list(product(*[list(range(len(i))) for i in values]))
        else:
            msg = "Could not determined range selection scheme from `concentrations`"
            raise ValueError(msg)

        self._validate_ranges(concentrations)

        N_states = len(self.c.states)
        # shape of [N_states, <values ligand 1>, <values ligand 2>, ..., <values ligand N>]
        res = np.zeros([N_states] + [len(i) for i in values])
        res_prob = np.zeros_like(res)
        res_std_error = np.zeros_like(res)
        res_covar = np.zeros([N_states, N_states] + [len(i) for i in values])
        res_deltas = np.zeros_like(res_covar)
        res_dGs = np.zeros_like(res_covar)

        pH = 0

        for g, c in zip(grid, indices):
            # using OrderedDict as a precaution mostly since
            # the shapes of the output can vary wildly
            conc = OrderedDict()
            coords = []

            for i, n in enumerate(names):
                if n == 'pH' or n.upper() == 'H+':
                    pH = g[i]
                    coords.append(i)
                    continue
                conc[n] = g[i]
                coords.append(i)

            self.c.concentrations = conc
            self.c.build_cycle(pH=pH)
            self.c.MLE(svd=svd)
            _filter1D = (slice(0, None), *c)
            _filter2D = (slice(0, None), slice(0, None), *c)
            res[_filter1D] = self.c.g_mle
            weights = np.exp(-res[_filter1D])
            Z = weights.sum()
            res_prob[_filter1D] = weights / Z
            res_std_error[_filter1D] = self.c.std_errors
            res_covar[_filter2D] = self.c.covariance_matrix
            res_deltas[_filter2D] = self.c.deltas
            res_dGs[_filter2D] = self.c.dGs

        coords = OrderedDict()
        coords['state'] = self.c.states.values[:, 0]
        coords['state_i'] = self.c.states.values[:, 0]
        coords['state_j'] = self.c.states.values[:, 0]

        for k, v in zip(names, values):
            coords[k] = v

        self.results = Dataset(
            data_vars=dict(
                free_energy=(
                    ['state', *names], res
                ),
                microstate_probs=(
                    ['state', *names], res_prob
                ),
                std_errors=(
                    ['state', *names], res_std_error
                ),
                covariance=(
                    ['state_i', 'state_j', *names], res_covar
                ),
                deltas=(
                    ['state_i', 'state_j', *names], res_deltas
                ),
                dGs=(
                    ['state_i', 'state_j', *names], res_dGs
                ),
            ),
            coords=coords,
        )

    def effective_energy_difference(self, category: str, state1: str, state2: str, **kwargs) -> float:
        self.c.g_mle = self.results.free_energy.sel(**kwargs).values
        return self.c.effective_energy_difference(category, state1, state2)

    def _validate_ranges(self, concentrations: dict) -> None:
        """Used to check the provided concentrations before executing the `run` method.

        Raises
        ------
        InvalidConcentrationError
        """

        for k, v in concentrations.items():
            for i in v:
                # since H+ is interpreted as a pH, we have to include all values
                if k in ('H+', 'pH'):
                    continue
                if i <= 0:
                    raise InvalidConcentrationError


class Multibind(object):
    """Core multibind object that constructs the potential graph and solves for thermodynamically
    consistent state free energies.

    Attributes
    ----------
    deltas : ndarray
        NxN numpy array with index i and j corresponding to the input free energy difference
        between states i and j. `np.nan` indicates an edge the was not defined by the user.
        This attribute becomes available after running the `MLE` method.
    covariance_matrix : ndarray
        NxN numpy array with index i and j corresponding to the covariance between the calculated
        state free energies. This attribute becomes available after running the `MLE` method.
    std_errors : ndarray
        1D numpy array where index i is the standard error of the calculated state free energies.
        This attribute becomes avaulable after running the `MLE` method.

    """

    def __init__(self, states_filename: str = None, graph_filename: str = None, comment_char: str = '#'):
        """
        Parameters
        ----------
        states_filename : str (optional)
            Path the CSV containing the states of the graph.
        graph_filename : str (optional)
            Path to the CSV containing the graph data for the network.
        comment_char : str, optional
            If there are comments in either files, you can specify the character. Defaults to '#'.

        """
        if states_filename:
            self.read_states(states_filename, comment=comment_char)
        else:
            self.states = None

        if graph_filename:
            self.read_graph(graph_filename, comment=comment_char)
        else:
            self.graph = None

        self.cycle = None
        self.concentrations = {}

    def build_cycle(self, pH: float = 5) -> None:
        """Constructs the cycle used for calculation

        Parameters
        ----------
        pH : float | int
            The pH to calculate binding free energies over.

        Returns
        -------
        None
        """

        # Determine if we have enough information to continue,
        # i.e. states information and graph information
        if self.states is None or self.graph is None:
            msg = "Need to specify both the state and graph \
            information. Try using `read_states` and `read_graph`."
            raise RuntimeError(msg)

        # Select all ligands that are not H+ and check if their concentrations
        # have been defined in the concentrations dictionary
        ligands = np.array(self.graph.ligand[(self.graph.ligand != "helm")
                                             & (self.graph.ligand != "H+")
                                             & (self.graph.ligand != "h+")])
        ligand_map = [x in self.concentrations.keys() for x in ligands]
        # if there are undefined ligand concentrations, raise an error and
        # warn the user
        if not all(ligand_map):
            missing_ligands = ligands[[not i for i in ligand_map]]
            msg = "Missing ligand concentrations for: {0}\n".format(" ".join(missing_ligands))
            msg += "Set them using the `concentrations` attribute"
            raise RuntimeError(msg)

        G = nx.DiGraph()
        # All states are defined in the states dataframe, use them for nodes
        G.add_nodes_from(self.states.name)

        # iterate through all connections
        for i in self.graph.index:
            # unpack for readability
            state1, state2, value, variance, ligand, standard_state = self.graph.iloc[i]

            # if we have protons, it must be a pKa
            if ligand.lower() == "h+":
                energy = protonation_free_energy(value, pH=pH)
                var = protonation_free_energy_standard_error(np.sqrt(variance))**2
            # using a direct helmholtz free energy
            elif ligand == "helm":
                energy = value
                var = variance
            # dealing with binding energies
            else:
                energy = binding_free_energy_general(value, concentration=self.concentrations[ligand])
                var = variance  # already in kT!
            # add a forward and direction energy

            G.add_edge(state1, state2, energy=energy, weight=var)
            G.add_edge(state2, state1, energy=-energy, weight=var)

        self.cycle = G

    def MLE(self, svd: bool = True) -> None:
        """Performs a maximum likelihood estimation on the current graph"""

        N = len(self.states.name)

        if svd:
            self._MLE_SVD(N)
        else:
            self._MLE_NR(N)

        self.std_errors = np.sqrt(np.diagonal(self.covariance_matrix))
        assert np.all(self.std_errors >= 0), "Standard errors for state free energies should be positive"

        self.g_mle = self.MLE_res - self.MLE_res[0]
        self.prob_mle = pd.DataFrame(np.exp(-self.g_mle) / np.sum(np.exp(-self.g_mle)), columns=["probability"])
        self.prob_mle["name"] = self.states.name

        self.deltas = np.zeros((N, N))
        self.dGs = np.zeros((N, N))

        self.deltas[:] = np.nan

        for r in self.graph.index:
            state1, state2, _, _, _, _ = self.graph.iloc[r]
            i = self.states[self.states.name == state1].index[0]
            j = self.states[self.states.name == state2].index[0]

            edge_attr = self.cycle.edges()[(state1, state2)]
            deltaij = edge_attr['energy']  # measured difference

            self.deltas[i, j] = deltaij
            self.deltas[j, i] = -deltaij

        for _i, dGi in enumerate(self.g_mle):
            for _j, dGj in enumerate(self.g_mle):
                self.dGs[_i, _j] = self.g_mle[_j] - self.g_mle[_i]

    def _MLE_SVD(self, N: int) -> None:
        B = np.zeros((N))
        A = np.zeros((N, N))

        for r in self.graph.index:
            state1, state2, _, _, _, _ = self.graph.iloc[r]
            i = self.states[self.states.name == state1].index[0]
            j = self.states[self.states.name == state2].index[0]

            edge_attr = self.cycle.edges()[(state1, state2)]
            deltaij = edge_attr['energy']  # measured difference
            varij = edge_attr['weight']  # measured variance

            B[i] += -deltaij / varij
            B[j] += deltaij / varij

            A[i, i] += 1 / varij
            A[j, j] += 1 / varij
            A[i, j] += -1 / varij
            A[j, i] += -1 / varij

        # since A is the Fisher information matrix, it's inverse is the covariance matrix of the
        # state free energies
        A_inv = pinv(A, hermitian=True)
        self.covariance_matrix = A_inv
        self.MLE_res = A_inv @ B

    def _MLE_NR(self, N: int) -> None:

        def grad_log_likelihood(g_t):
            """Returns the gradient of the log likelihood function.

            g_t : array of theoretical values for g
            """
            # state vector [g1, g2, g3, ... , gn-1, gn]
            state_vector = np.zeros(N)

            # factor that will be added to one node and subtracted from another
            def alphaij(gj, gi, deltaij, varij):
                return ((gj - gi) - deltaij) / varij

            # indices of state vector
            # Iterate over all connections
            for r in self.graph.index:
                state1, state2, value, variance, ligand, standard_state = self.graph.iloc[r]
                i = self.states[self.states.name == state1].index[0]
                j = self.states[self.states.name == state2].index[0]

                gj = g_t[j]
                gi = g_t[i]

                edge_attr = self.cycle.edges()[(state1, state2)]
                deltaij = edge_attr['energy']  # measured difference
                varij = edge_attr['weight']  # measured variance

                shared_alpha = alphaij(gj, gi, deltaij, varij)

                state_vector[i] += shared_alpha
                state_vector[j] -= shared_alpha

            return state_vector

        def kd(i, j):
            return int(i == j)

        def jacobian(g_t):
            # g_t here is not used deliberately as it is actually not needed except to avoid throwing an error
            J = np.zeros((N, N))
            for n in range(N):  # component of f
                for m in range(N):  # derivative with g_m
                    for k in self.graph.index:  # sum over ij
                        state1, state2, value, variance, ligand, standard_state = self.graph.iloc[k]

                        edge_attr = self.cycle.edges()[(state1, state2)]
                        varij = edge_attr['weight']  # measured variance

                        i = self.states[self.states.name == state1].index[0]
                        j = self.states[self.states.name == state2].index[0]
                        kdelta_factor = kd(n, j) * kd(m, i) - kd(n, j) * kd(m, j) - kd(n, i) * kd(m, i) + kd(n, i) * kd(m, j)
                        J[n, m] -= 1 / varij * kdelta_factor
            return J
        from scipy.optimize import root

        # use dijkstra_path to get the initial guess
        self.initial_guess = np.zeros(N)
        for i in range(1, N):
            edge_energies = nx.get_edge_attributes(self.cycle, 'energy')
            # edge_var = nx.get_edge_attributes(self.cycle, 'weight')
            path = nx.dijkstra_path(self.cycle, self.states.name[0], self.states.name[i])
            linked = [(path[j], path[j + 1]) for j, _ in enumerate(path[:-1])]
            self.initial_guess[i] = sum([edge_energies[x] for x in linked])

        FI = jacobian(None).T
        self.covariance_matrix = pinv(FI, hermitian=True)
        self.MLE_res = root(grad_log_likelihood, self.initial_guess, jac=jacobian).x

    def effective_energy_difference(self, macrostate_class: str, state1: str, state2: str) -> float:
        """Calculate the effective binding energy between two states.

        Parameters
        ----------
        macrostate_class : str
            Name of macrostate class (i.e. number of protons)
        state1 : str
            first, 'starting' state
        state2 : str
            second, 'destination' state

        Returns
        -------
        (float, float) : macroscopic free energy difference in kT and its associated error
        """

        macrostate_class = str(macrostate_class)

        microstates_1_indices = self.states[self.states[macrostate_class] == state1].index
        microstates_2_indices = self.states[self.states[macrostate_class] == state2].index

        energies_1 = self.g_mle[microstates_1_indices]
        energies_2 = self.g_mle[microstates_2_indices]

        std_err_1 = self.std_errors[microstates_1_indices]
        std_err_2 = self.std_errors[microstates_2_indices]

        z1 = np.exp(-energies_1)
        z2 = np.exp(-energies_2)

        std_err = np.sqrt((std_err_1**2 * z1**2).sum() / z1.sum()**2 + (std_err_2**2 * z2**2).sum() / z2.sum()**2)
        diff = np.log(np.sum(z1) / np.sum(z2))

        return diff, std_err

    def _parse(self, filename: str, comment: str = None) -> pd.DataFrame:
        """Helper function to quickly parse CSV into a DataFrame"""
        try:
            return pd.read_csv(filename, comment=comment)
        except pd.errors.EmptyDataError:
            mesg = f'Provided data file \'{filename}\' was empty'
            raise ValueError(mesg)
        except Exception as e:
            print(f'Could not parse file {filename}')
            raise e

    def read_states(self, filename: str, comment: str = None) -> pd.DataFrame:
        """Read in state information from a state CSV file.

        Parameters
        ----------
        filename : str
            Path to the state CSV file.

        Returns
        -------
        None
        """
        self.states = self._parse(filename, comment=comment)
        self.states['name'] = self.states['name'].astype('str')

    def read_graph(self, filename: str, comment: str = None):
        """Read in the graph information from a graph CSV file.

        Parameters
        ----------
        filename : str
            File path of the graph CSV

        Returns
        -------
        None
        """
        self.graph = self._parse(filename, comment=comment)
        self.graph['state1'] = self.graph['state1'].astype('str')
        self.graph['state2'] = self.graph['state2'].astype('str')
