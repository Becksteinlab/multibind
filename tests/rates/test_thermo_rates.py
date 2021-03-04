import pytest
from multibind.nonequilibrium import rate_matrix
import numpy as np
import os


class TestThermoRates(object):

    def test_matrix_regression(self):
        """Are results equal to the previous 'correct' results?"""
        c, matrix = rate_matrix("../examples/rates/inputs/inputs/rates.csv")
        regression_matrix = np.genfromtxt("../examples/outputs/outputs/test_rates.csv", delimiter=',')
        assert np.allclose(matrix, regression_matrix, atol=1e-5)