import pandas as pd
import numpy as np
import networkx as nx

k = 1.38064852e-23
avo = 6.023e23
R = k*avo / 1000.
T = 310.0

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
        self.concentations = {}

    def build_cycle(self, pH=5):
        """Constructs the cycle used for calculation"""

        # Determine if we have enough information to continue,
        # ie states information and graph information
        if not all([self.states, self.graph]):
            msg = "Need to specify both the state and graph
            information. Try using `read_states` and `read_graph`."
            raise ValueError(msg)

        # Select all ligands that are not H+ and check if their concentrations
        # have been defined in the concentations dictionary
        ligands = np.array(self.graph.ligand[(self.graph.ligand != "H+") & (self.graph.ligand != "h+")])
        ligand_map = [x in self.concentations.keys() for x in ligands]
        # if there are undefined ligand concentations, raise an error and
        # warn the user
        if not all(ligand_map):
            missing_ligands = ligands[[not i for i in ligand_map]]
            msg =  "Missing ligand concentations for: {0}\n".format(" ".join(missing_ligands))
            msg += "Set them using the `concentations` attribute"
            raise AttributeError(msg)

        G = nx.DiGraph()
        # All states are defined in the states dataframe, use them for nodes
        G.add_nodes_from(self.states.name)

        # iterate through all connections
        for i in self.graph.index():
            # unpack for readability
            state1, state2, value, variance, ligand, standard_state = self.graph.iloc[i]
            # if we have protons, it must be a pKa
            if ligand.lower() == "h+":
                energy = np.log(10)*(pH-value))
                var = np.log(10)**2 * variance
            # dealing with binding energies
            else:
                energy = value - np.log(self.concentrations[ligand]/standard_state)
                var = variance # already in kT!
            # add a forward and direction energy
            G.add_edge(state1, state2, energy=energy, weight=var)
            G.add_edge(state2, state1, energy=-energy, weight=var)

        self.cycle = G
            
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

    
