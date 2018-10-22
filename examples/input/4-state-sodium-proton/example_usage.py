import numpy as np
import multibind as mb
import matplotlib.pyplot as plt

def get_binding_energy(pH, Na):
    c = mb.Multibind(states_filename="states.csv", graph_filename="graph.csv")
    c.concentrations["Na+"] = Na
    c.build_cycle(pH=pH)
    c.MLE()
    return c.effective_energy_difference("bound","unbound","bound")

ph = np.linspace(1,14)

cs = [get_binding_energy(i,1) for i in ph]

plt.plot(ph,cs)
plt.show()
