# Installing multibind

## Install from PyPI

For most users, install the latest release with pip:

```bash
python -m pip install --upgrade pip
pip install multibind
```

## Requirements

- **Python** 3.10, 3.11, 3.12, 3.13, or 3.14 (see `pyproject.toml` for the exact supported range).
- **Runtime dependencies** (pulled in automatically): NumPy, pandas, SciPy, NetworkX, and xarray. NumPy 1.x is allowed on older Pythons; very new interpreters (for example 3.14) typically need a **NumPy 2.x** build with wheels for that version.

## Conda or Mamba

Create an environment with a suitable Python version, activate it, then use **pip** to install **multibind** as above. Avoid mixing multiple installers for the same package in one environment unless you know how they interact.


## For developers and contributors

The sections below cover working from a **Git clone**: versioning, Poetry, editable installs, tests, and documentation.

### Versioning and Git tags

Release versions on PyPI are produced with [versioningit](https://github.com/jwodder/versioningit) at **build** time from Git. Release tags must look like **`vMAJOR.MINOR.PATCH`** (for example `v0.2.0`); the leading `v` is removed in the published version string (see `pyproject.toml`).

Install from a **full** clone (`git clone`, not a tree without `.git`) when you build wheels or sdists locally so versioningit can read the repository.

### Install with Poetry

[Poetry](https://python-poetry.org/docs/#installation) resolves dependencies from `pyproject.toml` and is convenient for day-to-day development on a clone.

From the repository root:

```bash
poetry install --extras dev
```

**Why `--extras dev`?** The **`dev`** optional dependency group adds tools you need on a clone: **pytest**, **pytest-cov**, **Sphinx**, **IPython**, **versioningit**, and similar. Without it you cannot run the test suite or build the HTML docs in the usual way.

**Why is versioningit in `dev`?** The package version is **dynamic** (from Git tags via setuptools + versioningit when you run `pip install` / `python -m build`). Poetry’s **editable** install records a **`0.0.0`** placeholder in environment metadata. With **versioningit** installed, `multibind.__version__` can fall back to a real version derived from the repo (see `multibind/__init__.py`). Installs from PyPI always get the correct version from package metadata.

Runtime-only dependencies (no tests or docs tooling):

```bash
poetry install
```

Run the tests (paths assume the working directory is `tests/`):

```bash
cd tests
poetry run pytest -v
```

### Install with pip from a clone

With `.git` present so versioningit can run during the build:

```bash
python -m pip install --upgrade pip
pip install .
```

Editable install while you change the code (includes dev tools for tests/docs):

```bash
pip install -e ".[dev]"
```

### Install from GitHub

Replace the URL if you use a fork:

```bash
pip install "multibind @ git+https://github.com/BecksteinLab/multibind.git"
```

For a specific branch or tag, use the corresponding revision in the URL (see [pip VCS support](https://pip.pypa.io/en/stable/topics/vcs-support/)).

### Documentation

Online: [multibind.readthedocs.io](https://multibind.readthedocs.io/).

To build HTML locally, install dev dependencies (for example `poetry install --extras dev` or `pip install -e ".[dev]"`), then:

```bash
cd docs
sphinx-build -b html . _build/html
```
