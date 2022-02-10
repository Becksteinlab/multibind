import numpy.testing as npt
import multibind as mb
import numpy as np
from math import isclose

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
        self.c.MLE()
        # Found from maxima script
        # mb will output this order because of the state file order [g1, g2, g4, g3]
        dg_analytic = -8.73532566405552
        s1_indices = [0, 3]
        s2_indices = [1, 2]

        energies1 = self.c.g_mle[s1_indices]
        energies2 = self.c.g_mle[s2_indices]

        z1 = np.exp(-energies1)
        z2 = np.exp(-energies2)

        err_1 = self.c.std_errors[s1_indices]
        err_2 = self.c.std_errors[s2_indices]

        std_err_manual = np.sqrt(err_1[0]**2 * z1[0]**2 / z1.sum()**2
                                 + err_1[1]**2 * z1[1]**2 / z1.sum()**2
                                 + err_2[0]**2 * z2[0]**2 / z2.sum()**2
                                 + err_2[1]**2 * z2[1]**2 / z2.sum()**2
                                 )

        diff, std_err = self.c.effective_energy_difference("bound", "unbound", "bound")

        npt.assert_almost_equal(diff, dg_analytic, decimal=5)
        assert isclose(std_err, std_err_manual, rel_tol=1e-5)
