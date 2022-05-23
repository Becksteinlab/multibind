from importlib.metadata import version
from .multibind import Multibind, MultibindScanner, InvalidConcentrationError

__version__ = version("multibind")

del version
