from collections import OrderedDict
from multibind import InvalidConcentrationError
from multibind import Multibind
from multibind import MultibindScanner
import numpy as np
import numpy.testing as npt


class TestG_missing_coordinate(object):

    def setup(self):

        self.concentrations = OrderedDict(
            [('H+', [1, 2, 3, 4, 5, 6, 7]),
             ('Na+', [1, 0.250, 0.150, 0.100]),
             ]
        )

        self.scanner = MultibindScanner('../examples/input/4-state-diamond/states.csv',
                                        '../examples/input/4-state-diamond/graph.csv')

    def test_run(self):
        self.scanner.run(self.concentrations)
        g_analytic = np.array([4.237564495703078, 1.987494501321176,
                               0.0, -2.35510019160619])
        g_analytic = g_analytic - g_analytic[0]

        # since the graph doesn't have any Na+ dependence, we don't expect the answer to change
        npt.assert_almost_equal(self.scanner.results.free_energy.values[:, 4, 0], g_analytic)
        npt.assert_almost_equal(self.scanner.results.free_energy.values[:, 4, 1], g_analytic)
        npt.assert_almost_equal(self.scanner.results.free_energy.values[:, 4, 2], g_analytic)
        npt.assert_almost_equal(self.scanner.results.free_energy.values[:, 4, 3], g_analytic)

    def test_shape(self):
        """Does the algorithm give us the right number of results..."""
        self.scanner.run(self.concentrations)
        assert self.scanner.results.free_energy.values.shape == (4, 7, 4)


class TestG(object):

    def setup(self):

        self.concentrations = OrderedDict()

        self.concentrations['H+'] = np.linspace(0, 14, 5)
        self.concentrations['Na+'] = np.linspace(1, 0.001, 5)

        self.scanner_svd = MultibindScanner("../examples/input/4-state-sodium-proton/states.csv",
                                            "../examples/input/4-state-sodium-proton/graph.csv")

        self.scanner_NR = MultibindScanner("../examples/input/4-state-sodium-proton/states.csv",
                                           "../examples/input/4-state-sodium-proton/graph.csv")

        self.scanner_NR.run(self.concentrations, svd=False)
        self.scanner_svd.run(self.concentrations)

    def test_algorithm_equality(self):
        npt.assert_almost_equal(self.scanner_NR.results.free_energy.values,
                                self.scanner_svd.results.free_energy.values)
        npt.assert_almost_equal(self.scanner_NR.results.covariance.values,
                                self.scanner_svd.results.covariance.values)


def test_invalid_concentration():

    concentrations = OrderedDict()

    concentrations['H+'] = np.linspace(0, 14, 5)
    concentrations['Na+'] = np.linspace(-1, 0.1, 5)

    scanner = MultibindScanner("../examples/input/4-state-sodium-proton/states.csv",
                               "../examples/input/4-state-sodium-proton/graph.csv")

    try:
        scanner.run(concentrations)
    except InvalidConcentrationError:
        pass


def test_concentration_pairs():

    concentrations = [(1, 0.250), (2, 0.150), (3, 0.100)]

    scanner = MultibindScanner("../examples/input/4-state-sodium-proton/states.csv",
                               "../examples/input/4-state-sodium-proton/graph.csv")

    try:
        scanner.run(concentrations)
        assert False, "This should currently be failing."
    except NotImplementedError:
        pass


def test_missing_concentrations():

    scanner = MultibindScanner("../examples/input/4-state-sodium-proton/states.csv",
                               "../examples/input/4-state-sodium-proton/graph.csv")

    try:
        scanner.run(None)
        assert False, "Run function should throw an error when no concentrations are passed to it"
    except ValueError:
        pass


def test_effective_energy_difference():

    concentrations = {}
    concentrations["H+"] = [0, 7, 14]
    concentrations["Na+"] = [0.1, 0.5, 1]

    state_file = "../examples/input/4-state-sodium-proton/states.csv"
    graph_file = "../examples/input/4-state-sodium-proton/graph.csv"

    scanner = MultibindScanner(state_file, graph_file)
    c = Multibind(state_file, graph_file)

    scanner.run(concentrations)

    for n in concentrations["Na+"]:
        for h in concentrations["H+"]:
            c.concentrations["Na+"] = n
            c.build_cycle(pH=h)
            c.MLE()
            scanner_dG = np.array(scanner.effective_energy_difference("bound", "unbound", "bound", **{"Na+": n, "pH": h}))
            mb_dG = np.array(c.effective_energy_difference("bound", "unbound", "bound"))
            npt.assert_almost_equal(scanner_dG, mb_dG)
