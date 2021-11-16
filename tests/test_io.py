import pytest
import multibind as mb


def test_empty_state_file():
    try:
        mb.Multibind(states_filename='input/io_testing/empty.csv')
        assert False
    except ValueError:
        pass


def test_no_exist_file():
    try:
        mb.Multibind(states_filename='input/io_testing/no_exist.csv')
        assert False
    except FileNotFoundError:
        pass
