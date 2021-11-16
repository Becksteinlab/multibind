import pandas as pd
import numpy as np
import networkx as nx
from random import randint
from tqdm import tqdm


class MultibindDriver(object):

    """Class to quickly run multibind over a range of pH values"""

    def __init__(self, multibind):
        """
        Parameters
        ----------
        multibind : Multibind
            A multibind object with its states and graph attributed defined.

        """
        if not type(multibind.states) is None and not type(multibind.graph) is None:
            self.multibind = multibind
        else:
            raise ValueError(
                "Multibind driver must be passed a Multibind object that has states and a graph file loaded.")

    def create_tensor(self, pH_array):
        """Create a tensor containing state free energies across a pH range.

        Parameters
        ----------
        pH_array : iterable
            An iterable containing the pH values to calculate binding free energies for.

        Returns
        -------
        None

        """
        num_states = self.multibind.states.name.shape[0]
        self.tensor = np.zeros((num_states, num_states, len(pH_array)))

        for i, p in enumerate(pH_array):
            self.multibind.build_cycle(pH=p)
            self.multibind.MLE()
            for j in range(self.tensor.shape[1]):
                self.tensor[j, :, i] = self.multibind.g_mle - self.multibind.g_mle[j]


class Multibind(object):

    def __init__(self, states_filename=None, graph_filename=None):
        """
        Parameters
        ----------
        states_filename : str (optional)
            Path the CSV containing the states of the graph.
        graph_filename : str (optional)
            Path to the CSV containing the graph data for the network.
        """
        if states_filename:
            self.read_states(states_filename)
        else:
            self.states = None

        if graph_filename:
            self.read_graph(graph_filename)
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

    def MLE(self):
        """Performs a maximum likelihood estimation on the current cycle"""

        from scipy.optimize import root

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

        self.MLE_res = root(grad_log_likelihood, self.initial_guess, jac=jacobian)
        self.g_mle = self.MLE_res.x - self.MLE_res.x[0]
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
