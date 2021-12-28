import pandas as pd
import numpy as np
from numpy.linalg import pinv
import networkx as nx
from random import randint
from tqdm import tqdm
from xarray import Dataset
from itertools import product
from collections import OrderedDict


class InvalidConcentrationError(Exception):
    pass


class MultibindScanner(object):

    """Run multibind across a range of concentrations."""

    def __init__(self, statefile, graphfile, comment_char='#'):
        self.c = Multibind(states_filename=statefile, graph_filename=graphfile, comment_char=comment_char)

    def run(self, concentrations: dict, svd=True):
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
        res = np.zeros([N_states] + [len(i) for i in values])
        res_prob = np.zeros_like(res)

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
            _filter = (slice(0, None), *c)
            res[_filter] = self.c.g_mle
            weights = np.exp(-res[_filter])
            Z = weights.sum()
            res_prob[_filter] = weights / Z

        coords = OrderedDict()
        coords['state'] = self.c.states.values[:, 0]

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
            ),
            coords=coords,
        )

    def effective_energy_difference(self, category, state1, state2, **kwargs):
        self.c.g_mle = self.results.free_energy.sel(**kwargs).values
        return self.c.effective_energy_difference(category, state1, state2)

    def _validate_ranges(self, concentrations):
        """Used to check the provided concentrations before executing the `run` method.

        """

        for k, v in concentrations.items():
            for i in v:
                # since H+ is interpreted as a pH, we have to include all values
                if k in ('H+', 'pH'):
                    continue
                if i <= 0:
                    raise InvalidConcentrationError


class Multibind(object):

    def __init__(self, states_filename=None, graph_filename=None, comment_char='#'):
        """
        Parameters
        ----------
        states_filename : str (optional)
            Path the CSV containing the states of the graph.
        graph_filename : str (optional)
            Path to the CSV containing the graph data for the network.
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

    def build_cycle(self, pH=5):
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
        if type(self.states) is None or type(self.graph) is None:
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
                energy = np.log(10) * (pH - value)
                var = np.log(10) ** 2 * variance
            # using a direct helmholtz free energy
            elif ligand == "helm":
                energy = value
                var = variance
            # dealing with binding energies
            else:
                energy = value - np.log(self.concentrations[ligand] / standard_state)
                var = variance  # already in kT!
            # add a forward and direction energy

            G.add_edge(state1, state2, energy=energy, weight=var)
            G.add_edge(state2, state1, energy=-energy, weight=var)

        self.cycle = G

    def MLE(self, svd=True):
        """Performs a maximum likelihood estimation on the current graph"""

        N = len(self.states.name)

        def kd(i, j):
            return int(i == j)

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

        def jacobian(g_t):
            # g_t here is not used deliberately as it is actually not needed except to avoid throwing an error
            J = np.zeros((N, N))
            for n in range(N):  # component of f
                for m in range(N):  # derivative with g_m
                    for k in self.graph.index:  # sum over ij
                        state1, state2, value, variance, ligand, standard_state = self.graph.iloc[k]
                        i = self.states[self.states.name == state1].index[0]
                        j = self.states[self.states.name == state2].index[0]
                        kdelta_factor = kd(n, j) * kd(m, i) - kd(n, j) * kd(m, j) - kd(n, i) * kd(m, i) + kd(n, i) * kd(m, j)
                        J[n, m] += 1 / variance * kdelta_factor
            return J

        # use dijkstra_path to get the initial guess
        self.initial_guess = np.zeros(N)
        for i in range(1, N):
            edge_energies = nx.get_edge_attributes(self.cycle, 'energy')
            # edge_var = nx.get_edge_attributes(self.cycle, 'weight')
            path = nx.dijkstra_path(self.cycle, self.states.name[0], self.states.name[i])
            linked = [(path[j], path[j + 1]) for j, _ in enumerate(path[:-1])]
            self.initial_guess[i] = sum([edge_energies[x] for x in linked])

        if svd:
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

            A_inv = pinv(A, hermitian=True)
            self.MLE_res = A_inv @ B
        else:
            from scipy.optimize import root
            self.MLE_res = root(grad_log_likelihood, self.initial_guess, jac=jacobian).x

        self.g_mle = self.MLE_res - self.MLE_res[0]
        self.mle_linear_distortion = self.g_mle - (self.initial_guess - self.initial_guess[0])
        self.prob_mle = pd.DataFrame(np.exp(-self.g_mle) / np.sum(np.exp(-self.g_mle)), columns=["probability"])
        self.prob_mle["name"] = self.states.name
        return self.MLE_res

    def MLE_dist(self, N_steps=int(1e6), nt=1):
        """Run Monte-Carlo steps to assess quality of MLE results.

        Parameters
        ----------
        N_steps : int
            The number of Monte-Carlo steps to perform.

        Returns
        -------
        ndarray with the distribution of free energy values for the states.

        """

        def potential(g_t):
            potential = 0
            # factor that will be added to one node and subtracted from another
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

                potential += - 1 / (2 * varij) * ((gj - gi) - deltaij) ** 2

            return potential

        def accept(ns, cs):
            potential_ns = potential(ns)
            potential_cs = potential(cs)

            diff = potential_cs - potential_ns
            prob = min([1, np.exp(-50 * diff)])
            return np.random.random_sample() <= prob

        def compute(self, N_steps=N_steps):
            current_state = self.g_mle.copy()
            new_state = current_state.copy()
            step = 1
            accepted = 0
            rejected = 0
            Nstates = len(new_state)
            dist = np.zeros((N_steps - 1, Nstates))
            pbar = tqdm(total=N_steps, position=0)
            while step < N_steps:
                # select random state to mutate
                state = randint(0, Nstates - 1)
                # mutate state
                disp = np.random.normal(0, 0.01)
                new_state[state] = new_state[state] + disp
                # accept/reject change
                if accept(new_state, current_state):
                    current_state = new_state.copy()
                    dist[step - 1] = current_state[:]
                    pbar.update(1)
                    step += 1
                    accepted += 1
                else:
                    new_state = current_state.copy()
                    rejected += 1
            pbar.close()
            print("Accepted: ", accepted)
            print("Rejected: ", rejected)
            return dist

        return compute(self)

    def effective_energy_difference(self, macrostate_class, state1, state2):
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
        float : macroscopic free energy difference in kT
        """

        macrostate_class = str(macrostate_class)

        microstates_1_indices = self.states[self.states[macrostate_class] == state1].index
        microstates_2_indices = self.states[self.states[macrostate_class] == state2].index

        energies_1 = np.array([self.g_mle[i] for i in microstates_1_indices])
        energies_2 = np.array([self.g_mle[i] for i in microstates_2_indices])

        return np.log(np.sum(np.exp(-energies_1)) / np.sum(np.exp(-energies_2)))

    def _parse(self, filename, comment=None):
        """Helper function to quickly parse CSV into a DataFrame"""
        try:
            return pd.read_csv(filename, comment=comment)
        except pd.errors.EmptyDataError:
            mesg = f'Provided data file \'{filename}\' was empty'
            raise ValueError(mesg)
        except Exception as e:
            print(f'Could not parse file {filename}')
            raise e

    def read_states(self, filename, comment=None):
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

    def read_graph(self, filename, comment=None):
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
