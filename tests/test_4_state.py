import pytest
import numpy.testing as npt
import multibind as mb
import numpy as np

class TestG(object):

    def setup(self):
        self.c = mb.Multibind()
        self.c.read_graph("../examples/input/4-state-diamond/graph.csv",comment="#")
        self.c.read_states("../examples/input/4-state-diamond/states.csv")

    def test_length(self):
        """Does the algorithm give us the right number of results..."""
        self.c.build_cycle()
        self.c.MLE()
        assert len(self.c.g_mle) == 4

    def test_pH_5(self):
        self.c.build_cycle(pH=5)
        self.c.MLE()
        # Found from maxima script
        # mb will output this order because of the state file order [g1, g2, g4, g3]
        g_analytic = np.array([4.237564495703078, 1.987494501321176,
                               0.0, -2.35510019160619])
        g_analytic = g_analytic - g_analytic[0]
        
        npt.assert_almost_equal(self.c.g_mle,g_analytic)
        
    def test_microstate_probability(self):
        self.c.build_cycle(pH=5)
        self.c.MLE()
        g_analytic = np.array([4.237564495703078, 1.987494501321176,
                               0.0, -2.35510019160619])
        g_analytic = g_analytic - g_analytic[0]

        partition_function = np.sum(np.exp(-g_analytic))
        numerators = np.exp(-g_analytic)
        P = numerators/partition_function
        calculated_P = np.exp(-self.c.g_mle)/np.sum(np.exp(-self.c.g_mle))
        npt.assert_almost_equal(self.c.prob_mle.probability,P)
            
        

        

        
