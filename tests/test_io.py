import multibind as mb


def test_missing_graph():
    try:
        mb.Multibind(states_filename="../examples/input/4-state-diamond/states.csv").build_cycle()
        assert False
    except RuntimeError:
        pass


def test_missing_states():
    try:
        mb.Multibind(graph_filename="../examples/input/4-state-diamond/graph.csv").build_cycle()
        assert False
    except RuntimeError:
        pass


def test_missing_graph_and_states():
    try:
        mb.Multibind().build_cycle()
        assert False
    except RuntimeError:
        pass


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
