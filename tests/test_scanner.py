import numpy.testing as npt
from multibind.multibind import MultibindScanner
from collections import OrderedDict
import numpy as np
from multibind.multibind import InvalidConcentrationError


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
