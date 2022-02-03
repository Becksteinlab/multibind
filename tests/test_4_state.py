import numpy.testing as npt
import multibind as mb
import numpy as np
import pandas as pd
from multibind.chem import protonation_free_energy

class TestG(object):

    def setup(self):
        self.c = mb.Multibind()
        self.statefile = "../examples/input/4-state-diamond/states.csv"
        self.graphfile = "../examples/input/4-state-diamond/graph.csv"
        self.c.read_graph(self.graphfile, comment="#")
        self.c.read_states(self.statefile)

    def test_length(self):
        """Does the algorithm give us the right number of results..."""
        self.c.build_cycle()
        self.c.MLE()
        assert len(self.c.g_mle) == 4

    def test_pH_5_svd(self):
        self.c.build_cycle(pH=5)
        self.c.MLE(svd=True)
        # Found from maxima script
        # mb will output this order because of the state file order [g1, g2, g4, g3]
        g_analytic = np.array([4.237564495703078, 1.987494501321176,
                               0.0, -2.35510019160619])
        g_analytic = g_analytic - g_analytic[0]

        npt.assert_almost_equal(self.c.g_mle, g_analytic)

    def test_pH_5_NR(self):
        self.c.build_cycle(pH=5)
        self.c.MLE(svd=False)
        # Found from maxima script
        # mb will output this order because of the state file order [g1, g2, g4, g3]
        g_analytic = np.array([4.237564495703078, 1.987494501321176,
                               0.0, -2.35510019160619])
        g_analytic = g_analytic - g_analytic[0]

        npt.assert_almost_equal(self.c.g_mle, g_analytic)

    def test_microstate_probability(self):
        self.c.build_cycle(pH=5)
        self.c.MLE()
        g_analytic = np.array([4.237564495703078, 1.987494501321176,
                               0.0, -2.35510019160619])
        g_analytic = g_analytic - g_analytic[0]

        partition_function = np.sum(np.exp(-g_analytic))
        numerators = np.exp(-g_analytic)
        P = numerators / partition_function
        npt.assert_almost_equal(self.c.prob_mle.probability, P)

    def test_covariance_regression(self):
        """Test on previous results for the covariance matrix."""
        self.c.build_cycle()
        self.c.MLE()

        expected = np.array([[0.88016159, 0.51507475, -0.74761414, -0.6476222],
                             [0.51507475, 0.67087614, -0.6476222, -0.53832869],
                             [-0.74761414, -0.6476222, 0.88016159, 0.51507475],
                             [-0.6476222, -0.53832869, 0.51507475, 0.67087614]])

        npt.assert_almost_equal(self.c.covariance_matrix, expected, decimal=8)

    def test_equal_covariance(self):
        self.c.build_cycle()
        self.c.MLE()
        SVD_covar = self.c.covariance_matrix[:]

        self.c.build_cycle()
        self.c.MLE(svd=False)
        NR_covar = self.c.covariance_matrix[:]

        npt.assert_almost_equal(NR_covar, SVD_covar)

    def test_deltas(self):
        self.c.build_cycle(pH=0)
        self.c.MLE()

        states = pd.read_csv(self.statefile, comment='#')
        states['name'] = states['name'].astype('str')

        graph = pd.read_csv(self.graphfile, comment='#')
        graph['state1'] = graph['state1'].astype('str')
        graph['state2'] = graph['state2'].astype('str')

        state_indices = {j: i for i, j in enumerate(list(states['name'].values))}

        deltas = self.c.deltas

        for entry in range(graph.shape[0]):
            s1, s2, val, var, lig, std = graph.loc[entry]
            i1 = state_indices[s1]
            i2 = state_indices[s2]

            assert deltas[i1, i2] == protonation_free_energy(val, pH=0)
