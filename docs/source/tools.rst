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

- ``calculate_UU_correlators_V_slice.jl``: Use this to calculate the
  correlators (2-point functions) along U direction, in the
  double-null coordinates U,V. You first need to run a semiclassical
  simulations, either with or without bakcreaction, and save all the
  quantum modes for some timesteps. You create a grid for the
  coordinates U,V in post-porcessing, by setting ``NU``, ``NV``. You
  should first check ``correlators_UV_slice.ipynb`` (or
  ``correlators_UV_slice_subcritical.ipynb``) to identify for which
  V slice you want to calculate the U-U correlator.

- ``calculate_VV_correlators_U_slice.jl``: Similar to
  ``calculate_UU_correlators_V_slice.jl``, but for the V-V correlator,
  along a U slice.

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
  type for some timesteps of the simulation. Assumes a fixed ``rmax``.

- ``correlators_tr.ipynb``: Visualise the quantum, equal time,
  correlators (2-point functions in the r-domain). You need a
  semiclassical simulation, and to store ``data`` and
  ``correlators``. It assumes that both data types are saved on the
  same timesteps (if not, you need to make some modifications). It
  also assumes that you run with fixed ``rmax``.

- ``correlators_UV_slice.ipynb``: Visualisation for the U-U and V-V
  correlators, for a supercritical case (an apparent horizon is
  formed). You need to calculate these, using
  ``calculate_UU_correlators_V_slice.jl`` and
  ``calculate_VV_correlators_U_slice.jl``.

- ``correlators_UV_slice_subcritical.ipynb``: Similar to
  ``correlators_UV_slice.ipynb``, but for subcritical examples (no
  apparent horizon is formed).
  
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

- ``semiclassical_backreact_mink_consistency_check.ipynb``: Perform a
  consistency check for the backreaction implementation, against
  Minkowski. It calculates the L2 norms of the error between the
  numerical solution and Minkowski, and plots them. It requires some
  simulations with backreaction, with increasing values of ``kmax``
  and ``lmax``, and possibly ``mPV``. You need to store ``data`` and
  run with fixed ``rmax``. It is suggested to run with the same
  resolution for all simulations and save ``data`` at the same
  timesteps.

- ``semiclassical_evol_backreact_vs_no_backreact.ipynb``: Compare a
  simulation with and without backreaction. Run the same setup apart
  from the backreaction (on and off), and save ``data``,
  ``bilinears``, and ``quantum``. You could omit the latter, and
  appropriately comment out certain parts of the notebook.

- ``semiclassical_evol_standing_waves_mink_norms.jl``: Calculate norms
  for convergence test of the quantum modes, in Minkowski (the
  analytical solution is standing waves). You need to run the same
  simulation with three different resolutions, and save ``quantum`` at
  the same timesteps in all simulations. The script calculates the
  norm including all quantum modes, as well as a norm with only a
  specific subset. The latter needs to be tuned appropriately,
  depending on your simulation. It supports parallelization.

- ``semiclassical_evol_standing_waves_mink_plots.ipynb``: Plot the
  norms for the standing wave convergence test. First, you need to run
  ``semiclassical_evol_standing_waves_mink_norms.jl``.

- ``semiclassical_evol.ipynb``: Visualise various quantities from a
  single semiclassical simulation. Assumes that you saved ``data``,
  ``quantum``, and ``bilinears`` at the same timesteps, and that you
  run with fixed ``rmax``.

- ``semiclassical_setup_check_ID.ipynb``: Use this to check if your
  setup for a semiclassical simulation with backreaction is good
  (there are certain conditions that should be satisfied). You do not
  need a prior simulation.

- ``stress-energy_tensor_tr.ipynb``: Visualise various components of
  the stress-energy tensor, both classical and quantum
  contributions. You need a single semiclassical simulations (with or
  without backreaction) and to store ``data`` and ``bilinears``. It
  assumes that both data types are saved for the same timesteps (if
  not, you need to make some modifications). It also assumes that you
  run with fixed ``rmax``.

- ``stress-energy_tensor_UV.ipynb``: Visualise the energy (classical
  and semiclassical) in the outgoing null direction. You need a
  semiclassical simulation (with or without backreaction) where you
  store ``data`` and ``bilinears`` at the same timesteps. It assumes a
  fixed ``rmax``.

..
   - ``calculate_correlators_UV.jl``: Use this to calculate the
     correlators (2-point functions) in double-null coordinates U,V. You
     first need to run a semiclassical simulations, either with or
     without bakcreaction, and save all the quantum modes for some
     timesteps. You create a grid for the coordinates U,V in
     post-porcessing, by tuning ``NU``, ``NV``. It might be useful to
     check first ``correlators_UV.ipynb`` (maybe you need to comment out
     some lines in there), to see what is a good choice of ``NU`` and
     ``NV``.
..
   - ``correlators_UV.ipynb``: Visualise the quantum correlators across
     null directions (incoming or outgoing). The most relevant is the U-U
     correlators (across the outgoing null direction, keep V fixed), but
     you can also see the V-V correlator. After you perform a
     semiclassical simulation (with or without backreaction), and save
     ``data`` and ``quantum`` on the same timesteps, you need to run
     ``calculate_correlators_UV.jl``, for the same values of ``NU`` and
     ``NV``, as here. By default, the V-V correlators are commented out
     (also in the calculation).
