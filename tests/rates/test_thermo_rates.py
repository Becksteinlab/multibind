from multibind.nonequilibrium import rate_matrix
from multibind.nonequilibrium import project_rates, dG_standard_error, kij_standard_error, kji_standard_error
import numpy as np
import pytest
from math import isclose


class TestThermoRates(object):

    def test_matrix_regression(self):
        """Are results equal to the previous 'correct' results?"""
        c, matrix, SE = rate_matrix("../examples/rates/inputs/inputs/rates.csv")

        expected = np.array([[0.00000000e+00, 5.57242011e-02, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 2.79375886e+04],
                             [1.75947730e+10, 0.00000000e+00, 1.02000002e+10, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
                             [0.00000000e+00, 3.76999993e+07, 0.00000000e+00, 2.21842196e+05, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
                             [0.00000000e+00, 0.00000000e+00, 1.75947730e+10, 0.00000000e+00, 6.58803908e+04, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
                             [0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 4.86045874e+04, 0.00000000e+00, 1.75947730e+10, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
                             [0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 5.55550074e+05, 0.00000000e+00, 1.86999997e+07, 0.00000000e+00, 0.00000000e+00],
                             [0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 4.61000007e+09, 0.00000000e+00, 1.75947730e+10, 0.00000000e+00],
                             [0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.39547871e-01, 0.00000000e+00, 5.43754064e+04],
                             [1.15443517e+05, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 5.32970460e+04, 0.00000000e+00]])

        assert np.allclose(matrix, expected)


@pytest.mark.parametrize("kjib,kijb,sjib,sijb,target_dG,sj,si,kji,kij,sji,sij",
                         [(70, 40, np.sqrt(20), np.sqrt(20), -1, 1, 1, 21.30529001482594, 57.91378269735128, 34.74749438345482, 125.0178900064783)])
def test_projection(kjib, kijb, sjib, sijb, target_dG, sj, si, kji, kij, sji, sij):
    _k_ji, _k_ij = project_rates(kijb, kjib, sijb, sjib, target_dG)
    _sji = kji_standard_error(kijb, kjib, sijb, sjib, target_dG, sj, si)
    _sij = kij_standard_error(kijb, kjib, sijb, sjib, target_dG, sj, si)

    assert isclose(_k_ji, kji)
    assert isclose(_k_ij, kij)
    assert isclose(_sji, sji)
    assert isclose(_sij, sij)
