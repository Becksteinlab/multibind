# Installing multibind

## Requirements

- **Python** 3.10, 3.11, 3.12, 3.13, or 3.14 (see `pyproject.toml` for the exact supported range).
- **Dependencies** (installed automatically with the methods below): NumPy, pandas, SciPy, NetworkX, and xarray. NumPy 1.x is allowed on older Pythons; very new interpreters (for example 3.14) typically need a **NumPy 2.x** build that provides wheels for that version.

## Install with Poetry (recommended for development)

[Poetry](https://python-poetry.org/docs/#installation) manages the environment and resolves versions from `pyproject.toml`.

From the root of a clone of this repository:

```bash
poetry install
```

This installs **multibind** in editable mode plus **dev** tools (pytest, Sphinx, IPython, etc.).

To install only runtime dependencies (no dev extras):

```bash
poetry install --without dev
```

Run the test suite (paths in the tests assume the working directory is `tests/`):

```bash
cd tests
poetry run pytest -v
```

## Install with pip

### From a local checkout

Create and activate a virtual environment, then from the repository root:

```bash
python -m pip install --upgrade pip
pip install .
```

Editable install while you change the code:

```bash
pip install -e .
```

### From GitHub

Replace the URL if you use a fork:

```bash
pip install "multibind @ git+https://github.com/BecksteinLab/multibind.git"
```

For a specific branch or tag, use the corresponding revision in the URL (see [pip VCS support](https://pip.pypa.io/en/stable/topics/vcs-support/)).

### From PyPI

If a release is published on PyPI:

```bash
pip install multibind
```

## Conda or Mamba environments

Create an environment with a suitable Python version, activate it, then use **pip** inside that environment as above (from a clone or from Git/PyPI). Avoid mixing Poetry and conda installs of the same package in one environment unless you know how they interact.

## Documentation

To build the HTML docs locally, install dev dependencies (e.g. `poetry install`), then:

```bash
cd docs
sphinx-build -b html . _build/html
```

Online documentation: [multibind.readthedocs.io](https://multibind.readthedocs.io/).
