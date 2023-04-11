# Multibind


[![DOI](https://zenodo.org/badge/301552078.svg)](https://zenodo.org/badge/latestdoi/301552078) ![Tests](https://github.com/BecksteinLab/multibind/actions/workflows/tests.yml/badge.svg?branch=develop) [![codecov](https://codecov.io/gh/Becksteinlab/multibind/branch/develop/graph/badge.svg?token=7T3Z19P4W5)](https://codecov.io/gh/Becksteinlab/multibind) [![docs](https://readthedocs.org/projects/multibind/badge/?version=latest)](https://multibind.readthedocs.io/en/latest/?badge=latest)


Thermodynamic cycles which are defined by a set of individually determined free energy differences are not guaranteed to be thermodynamically consistent.
This criterion is only satisfied when all the differences in a closed loop vanish.
Multibind is a Python package that allows for the combination of these differences, along with their variances, into a set of data-informed and thermodynamically consistent state free energies.
Additionally, multibind supports cycles whose free energies are dependent on multiple ligand concentrations.

## State definitions

States are minimally defined by a name and should be added to a csv file.
Macrostate classes can be added as additional columns with the name of the class being the column header.

```text
name,n_protons,n_sodium
1,0,0
2,0,1
3,1,0
4,1,1
```

## Graph definition

States are not enough to define the graph, you'll need edges. 
These edges are the free energy differences between states and are characterized by a specific type of process.

Multibind allows for three process type, which are specified under the ligand column:

- "H+": proton binding which takes as its free energy value, the pKa.
- general ligand: binding of a general ligand. The standard state free energy in kT. By defining the concentration to build the cycle, this free energy becomes concentration dependent.
- "helm": undefined process which takes a free energy directly in kT.

In defining the edges, always treat the free energies as going from state 1 to state 2 for that process.
For 'H+', state 1 should be the deprotonated state and state 2 is the protonated state.
For a general ligand, state 1 is the unbound state and state 2 is the bound state.
Failure to do so will assign the negative of the desired free energy to that edge.

```text
# graph.csv
state1,state2,value,variance,ligand,standard_state
1,2,-10,0.1,Na+,1
2,3,5.697116,0.1,H+,1
4,3,-6,0.1,Na+,1
1,4,7.434294,0.1,H+,1
```

## Example code

The standard `Multibind` object is used to solve the graph at one point in concentration space.

```python
import multibind as mb

concentrations = {'Na': 0.150}

c = mb.Multibind()
c.read_graph("examples/input/4-state-diamond/graph.csv",comment="#")
c.read_states("examples/input/4-state-diamond/states.csv")
c.concentrations = concentrations
c.build_cycle(pH=7)
c.MLE(svd=True) 
# compute the effective energy difference between two macrostates
# if you are unsure of the values to pass to the function, look
# at either your state file or the contents of c.states
c.effective_energy_difference("N_protons",1,2)
```

The `MultibindScanner` can be used to scan across a bunch of concentrations.

```python
import multibind as mb

state_file = 'examples/input/4-state-diamond/states.csv'
graph_file = 'examples/input/4-state-diamond/graph.csv'

scanner = mb.MultibindScanner(state_file, graph_file, comment_char='#')

concentrations = {'H+': [7, 8, 9], 'Na+': [0.200, 0.150, 0.100]}
# This will MLE calculations across 9 points in concentration space
# the results will be in an xarray Dataset that can be accessed with
# scanner.results
scanner.run(concentrations, svd=True)
```

## Citation

When using **multibind** in published works, please cite the following preprint:

Kenney, Ian Michael, and Oliver Beckstein. Thermodynamically Consistent Determination of Free Energies and Rates in Kinetic Cycle Models. 2023. doi:10.1101/2023.04.08.536126.