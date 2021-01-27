# Multibind

Thermodynamic cycles which are defined by a set of individually determined free energy differences are not guaranteed to be thermodynamically consistent.
This criterion is only satisfied when all the differences in a closed loop vanish.
Multibind is a Python package that allows for the combination of these differences, along with their variances, into a set of data-informed and thermodynamically consistent state free energies.
Additionally, multibind supports cycles whose free energies are dependent on multiple ligand concentrations.

## State definitions

States are minimally defined by a name and should be added to a csv file.
Macrostates can be added as additional columns with the name of the category being the column header.

```text
name,bound
1,unbound
2,bound
3,bound
4,unbound
```

## Graph definition

States are not enough to define the graph, you'll need edges. 
These edges are the free energy differences between states and are characterized by a specific type of process.

Multibind allows for three process type, which are specified under the ligand column:

- "H+": proton binding which takes as its free energy value, the pKa.
- "helm": undefined process which takes a free energy directly in kT.
- general ligand: binding of a general ligand. The standard state free energy in kT. By defining the concentration to build the cycle, this free energy becomes concentration dependent.

In defining the edges, always treat the free energies as going from state 1 to state 2.

```text
state1,state2,value,variance,ligand,standard_state
1,2,-10,0.1,Na+,1
2,3,5.697116,0.1,H+,1
4,3,-6,0.1,Na+,1
1,4,7.434294,0.1,H+,1
```

## Example code

```python
import multibind as mb

c = mb.Multibind()
c.read_graph("examples/input/4-state-diamond/graph.csv",comment="#")
c.read_states("examples/input/4-state-diamond/states.csv")
c.build_cycle(pH=7)
c.MLE() 
# compute the effective energy difference between two macrostates
# if you are unsure of the values to pass to the function, look
# at either your state file or the contents of c.states
c.effective_energy_difference("Nprotons",1,2)
```