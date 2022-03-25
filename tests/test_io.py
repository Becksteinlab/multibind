import multibind as mb
from pathlib import Path
import tempfile


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


def test_successful_read_modify_write_read():
    with tempfile.NamedTemporaryFile() as F:
        states = Path() / ".." / "examples" / "input" / "4-state-diamond" / "states.csv"
        graph = Path() / ".." / "examples" / "input" / "4-state-diamond" / "graph.csv"
        c1 = mb.Multibind(states_filename=states, graph_filename=graph)
        c1.graph.loc[1, ["value"]] = 8.0
        c1.write_graph(F.name)

        c2 = mb.Multibind(states_filename=states, graph_filename=F.name)

        assert c1.graph.loc[1, ["value"]].value == c2.graph.loc[1, ["value"]].value
