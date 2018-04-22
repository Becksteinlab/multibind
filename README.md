Example usage:

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