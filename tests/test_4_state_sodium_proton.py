import pytest
import numpy.testing as npt
import multibind as mb
import numpy as np

class TestG(object):

    def setup(self):
        self.c = mb.Multibind()
        self.c.read_graph("../examples/input/4-state-sodium-proton/graph.csv")
        self.c.read_states("../examples/input/4-state-sodium-proton/states.csv")
        self.c.concentrations["Na+"] = 1

    def test_length(self):
        """Does the algorithm give us the right number of results..."""
        self.c.build_cycle()
        self.c.MLE()
        assert len(self.c.g_mle) == 4

    def test_pH_7(self):
        self.c.build_cycle(pH=7)
        print(self.c.graph)
        self.c.MLE()
        # Found from maxima script
        # mb will output this order because of the state file order [g1, g2, g4, g3]
        dg_analytic = -8.73532566405552
        npt.assert_almost_equal(self.c.effective_energy_difference("bound","unbound","bound"),dg_analytic, decimal=5)
        

        
