from setuptools import setup, find_packages

setup(name="multibind",
      version="0.1.0",
      description="A library for calcuating effective free energies of systems with multiple protonation states.",
      author="Ian Kenney",
      author_email="ian.kenney@asu.edu",
      packages=find_packages(),
      install_requires=['numpy', 'networkx', 'pandas', 'matplotlib','scipy'],
      tests_require=['pytest']
)
