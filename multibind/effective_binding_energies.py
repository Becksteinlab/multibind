import itertools

import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

k = 1.38064852e-23
avo = 6.023e23
R = k*avo / 1000.
T = 310.0

# THIS IS PROBLEM SPECIFIC AND SHOULD BE MOVED EVENTUALLY
residues = {0: "asp156", 1: "asp157", 2: "lys305"}

def _links(x):
    """Helper function to link states with one another."""
    return [(x[i],x[i+1]) for i,_ in enumerate(x[:-1])]

def binding_energy(DG, Na, c0=1):
    return np.float128(DG - (R*T) * np.log(Na/c0))

def protonation_energy(pka, pH):
    return np.float128((R*T) * np.log(10)*(pH-pka))

def deprotonation_energy(pka, pH):
    return  -protonation_energy(pka, pH)

def get_path_energy(cycle, s1, s2):
    edge_energies = nx.get_edge_attributes(cycle, 'energy')
    edge_var = nx.get_edge_attributes(cycle, 'weight')
    path = nx.dijkstra_path(cycle, s1, s2)
    linked = _links(path)
    res_sum = sum([edge_energies[x] for x in linked])
    res_var = sum([edge_var[x] for x in linked])
    return res_sum, res_var, path

def energy_matrix(cycle):
    """Create the energy "matrix" G[i, j] as a DataFrame
    
    The free energy for going from state i to j is
    
         G.loc[i][j]
    """
    nodes = sorted(cycle.nodes())
    result = np.zeros((len(nodes), len(nodes)))
    result_var = np.zeros_like(result)
    for i in range(1, len(nodes)):
        for j in range(0, i):
            energy, variance, _ = get_path_energy(cycle, nodes[i], nodes[j])
            result[i,j] = energy
            result[j,i] = -energy
            result_var[i,j] = result_var[j,i] = variance       
    
    labels = [cycle.nodes.data('label')[node] for node in nodes]
    energies = pd.DataFrame(result, index=labels, columns=labels, dtype=np.float128)
    var = pd.DataFrame(result_var, index=labels, columns=labels, dtype=np.float128)
    
    return energies, var, labels

def general_binding(energy_matrix, var_matrix, reference_state='uS1'):
    bound_states = [label for label in energy_matrix.index if label.startswith('b')] 
    unbound_states = [label for label in energy_matrix.index if label.startswith('u')]
    
    beta = 1/(R*T)
    G = energy_matrix.loc[reference_state]
    E = var_matrix.loc[reference_state]
    Z = np.sum(np.exp(-beta * G))
    B = np.sum(np.exp(-beta * G.drop(unbound_states)))
    bound_factors = np.exp(-2 * beta * G.drop(unbound_states))/(B**2)
    bound_var = E.drop(unbound_states)
    bound_variance_cont = np.sum(bound_factors * bound_var)
    unbound_factors = np.exp(-2 * beta * G.drop(bound_states))/(Z-B)**2
    unbound_var = E.drop(bound_states)
    unbound_variance_cont = np.sum(unbound_factors * unbound_var)
    total_variance = unbound_variance_cont + bound_variance_cont
    return np.log(Z/B - 1), total_variance

def make_states(residues):
    assert len(residues) == 3  # should make 
    state_vectors = list(itertools.product(*(len(residues) * [[0, 1]])))
    labels = "S4 S2 8 S1 6 7 5 S3".split()
    states = {}
    for b_label in ("b", "u"):
        bstates = {b_label + label:(b_label, state) for 
               label, state in zip(labels, state_vectors)}
        states.update(bstates)
    return states

def protonation_DG(pkas, resname, pH):
    return [protonation_energy(pkas.loc[resname]['mean'], pH),
            (R*T)**2 * np.log(10)**2 * pkas.loc[resname]['var']]

def make_graph(pkas, binding, conformation="IF", pH=7, c_na=1.0,
               bound_states=(('b', (0,0,1)), ('b', (0,0,0)), # known: bS2, bS4
                             ('b', (1,0,0))),                # putative
                             # other states such as S3 (1,1,1) or (0,1,1) no binding in eqMD 
               residues=residues):
    states = make_states(residues)
    
    pka = pkas[conformation]

    G = nx.DiGraph()
    G.add_nodes_from((state, {"label": label}) for label, state in states.items())

    # add edges for protonation/deprotonation
    nodes = list(G.nodes.keys())
    for node in nodes:
        b_label, state = node
        for resindex, value in enumerate(state):
            newstate = np.array(state)
            if value == 0:
                # add a proton
                newstate[resindex] += 1
                factor = 1
            elif value == 1:
                # remove a proton
                newstate[resindex] -= 1
                factor = -1
            else:
                raise ValueError(("WTF: for state = {}, value = {} "
                                 "should not happen").format(state, value))
            newnode = (b_label, tuple(newstate))
            DG, varDG = protonation_DG(pka[b_label], residues[resindex], pH)
            DG *= factor  # -1 for deprotonation
            G.add_edge(node, newnode, energy=DG, weight=varDG)
            
    # remove bound nodes which we have not observed
    non_existant = [n for n in G if n[0] == "b" and n not in bound_states]
    G.remove_nodes_from(non_existant)

    # add edges for Na+ binding between u and b nodes
    binding_edges = [("uS2", "bS2"), ("uS4", "bS4")]
    # select data for conformation ....
    for b in binding_edges:
        label = (b[0][1:] + conformation).lower()  # eg uS2 --> s2if
        energy = binding.loc[label]['mean']
        var = binding.loc[label]['var']
        DGbind = binding_energy(energy, c_na)
        G.add_edge(states[b[0]], states[b[1]], 
                   weight=var, energy=DGbind)
        G.add_edge(states[b[1]], states[b[0]], 
                   weight=var, energy=-DGbind)
    return G, states

def calculate(pkas, binding, conformation="IF", pH=7, c_na=1.0, residues=residues):

    cycle, states = make_graph(pkas, binding, conformation=conformation,
                              pH=pH, c_na=c_na, residues=residues)
    matrix_results = energy_matrix(cycle)
    return matrix_results, cycle
