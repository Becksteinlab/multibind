import pandas as pd
import numpy as np
import networkx as nx

class MultibindDriver(object):

    def __init__(self, multibind):
        if not type(multibind.states) is None and not type(multibind.graph) is None:
            self.multibind = multibind
        else:
            raise ValueError("Multibind driver must be passed a Multibind object that has states and a graph file loaded.")

    def create_tensor(self, pH_array):
        num_states = self.multibind.states.name.shape[0]
        self.tensor = np.zeros((num_states, num_states, len(pH_array)))
        
        for i,p in enumerate(pH_array):
            self.multibind.build_cycle(pH=p)
            self.multibind.MLE()
            for j in range(self.tensor.shape[1]):
                self.tensor[j,:,i] = self.multibind.g_mle - self.multibind.g_mle[j]

            

class Multibind(object):

    def __init__(self, states_filename=None, graph_filename=None):
        # If states are specified in a CSV, may as well fill them in
        # here. The same goes for the graph information
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
        """Constructs the cycle used for calculation"""

        # Determine if we have enough information to continue,
        # ie states information and graph information
        if type(self.states) == None or type(self.graph) == None:
            msg = "Need to specify both the state and graph \
            information. Try using `read_states` and `read_graph`."
            raise RuntimeError(msg)

        # Select all ligands that are not H+ and check if their concentrations
        # have been defined in the concentrations dictionary
        ligands = np.array(self.graph.ligand[(self.graph.ligand != "H+") & (self.graph.ligand != "h+")])
        ligand_map = [x in self.concentrations.keys() for x in ligands]
        # if there are undefined ligand concentrations, raise an error and
        # warn the user
        if not all(ligand_map):
            missing_ligands = ligands[[not i for i in ligand_map]]
            msg =  "Missing ligand concentrations for: {0}\n".format(" ".join(missing_ligands))
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
                energy = np.log(10)*(pH-value)
                var = np.log(10)**2 * variance
            # dealing with binding energies
            else:
                energy = value - np.log(self.concentrations[ligand]/standard_state)
                var = variance # already in kT!
            # add a forward and direction energy

            G.add_edge(state1, state2, energy=energy, weight=var)
            G.add_edge(state2, state1, energy=-energy, weight=var)

        self.cycle = G

    def MLE(self):
        """Performs a maximum likelihood estimation on the current cycle"""

        from scipy.optimize import root

        N = len(self.states.name)
        kd = lambda i,j : int(i==j)
        
        def grad_log_likelihood(g_t):
            """Returns the gradient of the log likelihood function.
            
            g_t : array of theoretical values for g
            """
            # state vector [g1, g2, g3, ... , gn-1, gn]
            state_vector = np.zeros(N)

            # factor that will be added to one node and subtracted from another
            alphaij = lambda gj, gi, deltaij, varij: ((gj - gi) - deltaij) / varij

            # indices of state vector
            # Iterate over all connections
            for r in self.graph.index:
                state1, state2, value, variance, ligand, standard_state = self.graph.iloc[r]
                i = self.states[self.states.name == state1].index[0]
                j = self.states[self.states.name == state2].index[0]
                
                gj = g_t[j]
                gi = g_t[i]
                
                edge_attr = self.cycle.edges()[(state1,state2)]
                deltaij = edge_attr['energy'] # measured difference
                varij = edge_attr['weight'] # measured variance
                
                shared_alpha = alphaij(gj, gi, deltaij, varij)
                
                state_vector[i] += shared_alpha
                state_vector[j] -= shared_alpha
                
            return state_vector

        def jacobian(g_t):
            J = np.zeros((N,N))
            for n in range(N): # component of f
                for m in range(N): # derivative with g_m
                    for k in self.graph.index: # sum over ij
                        state1, state2, value, variance, ligand, standard_state = self.graph.iloc[k]
                        i = self.states[self.states.name == state1].index[0]
                        j = self.states[self.states.name == state2].index[0]
                        kdelta_factor = kd(n,j)*kd(m,i) - kd(n,j)*kd(m,j) - kd(n,i)*kd(m,i) + kd(n,i)*kd(m,j)
                        J[n,m] += 1/variance * kdelta_factor
            return J
        
        # use dijkstra_path to get the initial guess
        initial_guess = np.zeros(N)
        for i in range(1,N):
            edge_energies = nx.get_edge_attributes(self.cycle, 'energy')
            edge_var = nx.get_edge_attributes(self.cycle, 'weight')
            path = nx.dijkstra_path(self.cycle, self.states.name[0], self.states.name[i])
            linked = [(path[j],path[j+1]) for j,_ in enumerate(path[:-1])]
            initial_guess[i] = sum([edge_energies[x] for x in linked])

        self.MLE_res = root(grad_log_likelihood, initial_guess, jac=jacobian)
        self.g_mle = self.MLE_res.x - self.MLE_res.x[0]
        self.prob_mle = pd.DataFrame(np.exp(-self.g_mle)/np.sum(np.exp(-self.g_mle)),columns=["probability"])
        self.prob_mle["name"] = self.states.name
        return self.MLE_res

    def effective_energy_difference(self, macrostate_class, state1, state2):
        """Calculate the effective binding energy between two states.
        
        Parameters
        ==========
        macrostate_class : name of macrostate class (i.e. number of protons)
        state1 : first, 'starting' state
        state2 : second, 'destination' state

        Returns
        =======
        float : binding free energy in kT
        """

        macrostate_class = str(macrostate_class)
        
        microstates_1_indices = self.states[self.states[macrostate_class] == state1].index
        microstates_2_indices = self.states[self.states[macrostate_class] == state2].index

        energies_1 = np.array([self.g_mle[i] for i in microstates_1_indices])
        energies_2 = np.array([self.g_mle[i] for i in microstates_2_indices])

        return np.log(np.sum(np.exp(-energies_1))/np.sum(np.exp(-energies_2)))

    def _parse(self, filename, comment=None):
        """Helper function to quickly parse CSV into a DataFrame"""
        try:
            return pd.read_csv(filename, comment=comment)
        except Exception as e:
            raise e("Could not parse file %" % filename)
            
    def read_states(self, filename, comment=None):
        """Read in state information from a state CSV file

        Parameters
        ==========
        filename : string with the file path
        """
        self.states = self._parse(filename, comment=comment)
        self.states['name'] = self.states['name'].astype('str')
        
    def read_graph(self, filename, comment=None):
        """Read in the graph information from a graph CSV file

        Parameters
        ==========
        filename : string with the file path

        Returns
        =======
        DataFrame with graph information (accessible using `graph`
        attribute)
        """
        self.graph = self._parse(filename, comment=comment)
        self.graph['state1'] = self.graph['state1'].astype('str')
        self.graph['state2'] = self.graph['state2'].astype('str')
