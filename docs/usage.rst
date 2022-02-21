Usage
=====

Installation
------------

To use multibind, first install it using pip:

.. code-block:: console

   pip install multibind


Defining a graph
----------------

Graphs can be fully defined with two structures: 1) a state file containing state names and how those states group into macrostates and 2) a graph file defining how the states are connected and by which processes characterize the connections.
As an example consider a system that has two titratable sites, which we refer to as left and right.
Each site can have a protein bound or unbound.
This gives 4 possible states, whos labels can be defined by their binary representations for site binding:

    - 00 (left deprotonated, right deprotonated)
    - 01 (left deprotonated, right protonated)
    - 10 (left protonated, right deprotonated)
    - 11 (left protonated, right protonated)

With just this information, the state file can be defined as follows:

.. code-block::

   state,N
   00,0
   01,1
   10,1
   11,2

where we can include arbitrary columns defining potentially relevant macrostates.
In this case, we know how many protons are currently bound to the molecule, so we included the column `N` to signify the three macrostates of number of protons bound (0, 1, and 2).
Defining these now will make calculations such as macroscopic pKas easier after we've built the model.
Including states early on does not influence how the model is built, so add as many physically meaningful categories as you'd like.
Now, to define the graph file, the only things we need to know are: 1) which connections do we have information about and 2) what process characterizes the connection.
Here, we are only dealing with protonation and deprotonation events.
Let's suppose that we have calculated site specific pKas (`pX`) (along with error estimates in the form of variances, `vX`) through simulation for the following connections:

    - 00 → 01 (p1, v1)
    - 00 → 10 (p2, v2)
    - 01 → 11 (p3, v3)
    - 10 → 11 (p4, v4)

The graph file is then defined as:

.. code-block::

   state1,state2,value,variance,ligand,standard_state
   00,01,p1,v1,H+,1
   00,10,p2,v2,H+,1
   01,11,p3,v3,H+,1
   10,11,p4,v4,H+,1

`state1` and `state2` are the starting and ending states, respectively.
If a ligand other than `helm` is specified, then `state2` must always be the state that gained the binding ligand.
The `value` column, in general, is the free energy difference between the states at standard state concentration (1 M).
In the case where the ligand is specified as `H+`, the `value` column then takes a pKa since this is typically what is reported.
The `variance` column is simply the variance of the measured `value`.

.. list-table::
   :header-rows: 1

   * - Ligand
     - Description
   * - H+
     - The `value` should be the pKa describing the process. The `variance` is the variance in the measured `value`. The input for `standard_state` does not matter.
   * - helm
     - The `value` should be the fixed free energy between states (F\ :sub:`2` - F\ :sub:`1`). The `variance` is the variance of this value. The `standard_state` does not apply to this type of entry (typically left as 1 M).
   * - `*`
     - Anything other than the two previous rows is treated as an arbitrary ligand species. The `value` is the standard state free energy difference between the states. The `variance` is the variance of that `value`. The `standard_state` defines where in concentration space the `value` is valid. Typically, `standard_state` is just at 1 M.

Note that if instead of a proton, the right binding site bound to a different ligand, say a sodium ion, then the state file could have been defined as


.. code-block::

   state,N_prot,N_sod
   00,0,0
   01,0,1
   10,1,0
   11,1,1

and the graph files could have been defined as

.. code-block::

   state1,state2,value,variance,ligand,standard_state
   00,01,F1,v1,Na+,1
   00,10,p2,v2,H+,1
   01,11,p3,v3,H+,1
   10,11,F4,v4,Na+,1

where p1 (pKa 1) became F1 (free energy difference 1) and p4 (pKa 4) became F4 (free energy difference 4).
The `v1` and `v4` values now refer to the variance in the standard state sodium binding free energy difference for connections 1 and 4, respectively.
The ligand label also changed to reflect that this connection is now a function of the sodium concentration rather than the pH.

Solving for a graph
-------------------
Once the state and graph files have been defined, we are ready to build the cycle and solve for the free energies.
For a predefined graph, we will use the two site, four state proton binding model from above with a state and graph file defined as:

.. code-block::

   # states.csv
   name,Nprotons
   1,0
   2,1
   4,1
   3,2

.. code-block::

   # graph.csv
   state1,state2,value,variance,ligand,standard_state
   1,2,6,0.1,H+,1
   2,3,7,0.5,H+,1
   4,3,6,0.1,H+,1
   1,4,5.7,5,H+,1

The system can then be loaded and solved using the `Multibind` object.

.. code-block:: python

   >>> import multibind as mb
   >>> c = mb.Multibind(state_file="states.csv", graph_file="graph.csv", comment='#')
   >>> c.build_cycle(pH=8)
   >>> c.MLE()
   >>> c.g_mle # numpy array containing the state free energies at pH 8

The entries of the `g_mle` attribute correspond to the order of the states as defined in `states.csv`.

Macroscopic free energies can be calculated.
For instance, in order to calculate the macroscopic binding free energy between a zero bound proton state and a single bound proton state, you can use

.. code-block:: python

   >>> c.effective_energy_difference("Nprotons", 0, 1)

This binding free energy might a quantity of interest over a range of pH values.
Rather than iterating over a range of values, the `MultibindScanner` conveniently handles this for you and stores standard results in an `xarray.DataArray`.
For the above example over the pH range 0 to 14, we can write:

.. code-block:: python

   >>> from multibind.multibind import MultibindScanner
   >>> import numpy as np
   >>> scanner = MultibindScanner("states.csv", "graph.csv", comment_char='#')
   >>> concentrations = {"H+": np.linspace(0, 14)}

   >>> scanner.run(concentrations)
   >>> scanner.results.pH # pH values scanned over
   >>> scanner.results.microstate_probs # miscrostate probabilities as a function of pH
   >>> scanner.results.free_energy # state free energies as a function of pH
