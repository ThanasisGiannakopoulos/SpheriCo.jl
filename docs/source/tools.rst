.. _Tools:

Tools
=======

There are various tools under ``SpheriCo.jl/tools/`` to help you
analyse your simulations. All of them are written in Julia and most of
them are interactive `Jupyter notebooks <https://jupyter.org/>`_. Some
are Julia scripts that you can run e.g. in your terminal, like the
:ref:`examples <examples>`. Instructions on how to use each script are
provided within it. Below you can see briefly what each tool does and
how it interacts with other tools:

- ``bilinears_evol.ipynb``: Analyse the (potential) effect of the
  quantum modes on the classical geometry, both with or without
  (potential) backreaction. First, you need to run a semiclassical
  simulation and save the bilinears and classical variables. The
  script also loads the quantum variables, but are not used (so saving
  them is optional for this script, though you need to comment the
  relevant lines out).

- ``bilinears_ID.ipynb``: See the effect of the quantum modes to the
  matter content at the right-hand-side of the equations for K, KB
  (extrinsic curvature). You can tune the quantum numbers ``kmax``,
  ``lmax``, the Pauli-Villars mass ``mPV``, the convention for Planck
  mass, and more, and see how close the **initial** stress-energy
  tensor is to Minkowski, when the specific backreaction is
  considered. It might be useful to check this before running a
  semiclassical example with backreaction.

- ``calculate_UV_correlators.jl``: Use this to calculate the
  correlators (2-point functions) in double-null coordinates U,V. You
  first need to run a semiclassical simulations, either with or
  without bakcreaction, and save all the quantum modes for some
  timesteps. You create a grid for the coordinates U,V in
  post-porcessing, by tuning ``NU``, ``NV``. It might be useful to
  check first ``UV_correlators.ipynb`` (maybe you need to comment out
  some lines in there), to see what is a good choice of ``NU`` and
  ``NV``.

- ``check_Hamiltonian_momentum_reduction_constraints.ipynb``: Check
  the violation of the classical Hamiltonian, momentum and reduction
  constraints. You can compare this violation between different
  simulations (e.g. with and without constraint damping). It considers
  only the classical geometry, but it can also be used for a
  semiclassical simulation, as long as the classical data are stored.

- ``check_reduction_constraints_convergence.ipynb``: Inspect the
  pointwise convergence of the classical reduction variables. First,
  you need to run the same simulation in three different resolutions.

- ``classical_norm_convergence.ipynb``: Check for convergence in the
  L2 norm of the classical state vector, for strong and weak data. You
  need to simulate the same setup in three different resolutions
  (double ``Nr`` every time) and save ``data``.

- ``classical_norm_convergence_noisy.ipynb``: Perform norm convergence
  tests with noisy data. You need at least three different resolutions
  for the same setup, and to store ``data``.

- ``classical_pointwise_convergence.jl``: Test the pointwise
  convergence of various classical functions (in r). Assumes three
  simulations of the same setup, with the only difference the
  resolution (``Nr``). Also assumes that every time resolution is
  increased, ``Nr`` is doubled (increase ``D`` in the examples by
  one). You need to store ``data``. You also need to modify the
  script, and give the directory where your data are stored.

- ``classical_scalar_evolution.ipynb``: Analyse various classical
  quantities in a single simulation. You need to store the ``data``
  type for some timesteps of the simulation.

- ``correlators.ipynb``:
  
- ``criticality_echos_convergence.ipynb``: Examine the pointwise
  convergence at r=0 of the scalar field and lapse function, as well
  as the behaviour of the Ricci scalar and other quantities, for near
  critical simulations. It assumes that you run for at least three
  different resolutions and store the ``r0data`` type. It is suggested
  to run with infalling rmax, but fixed rmax is also an option
  (possibly not good results though).

- ``criticality_echos_infalling_rmax_evol.ipynb``: Analyse a single
  near critical simulation. You can plot the behaviour of the scalar
  field, lapse, Ricci scalar, and more, at r=0, as well as create
  movies of these quantities for the whole radial domain, in time. You
  need to save the data types ``r0data`` and ``data``,
  respectively. It works as is for infalling ``rmax``, but with
  minimal changes for fixed ``rmax`` as well.

- ``evol_quantum_backreact_vs_no_backreact.ipynb``:

- ``evol_quantum.ipynb``:

- ``evol_quantum_standing_waves_mink_norms.jl``:

- ``evol_quantum_standing_waves_mink_plots.ipynb``:

- ``mink_backreact_consistency_check.ipynb``:

- ``quantum_setup_check_ID.ipynb``:

- ``UV_correlators.ipynb``:
