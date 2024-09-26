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

- ``bilinears_evol.ipynb``: analyse the (potential) effect of the
  quantum modes on the classical geometry, both with or without
  (potential) backreaction. First, you need to run a semiclassical
  simulation and save the bilinears and classical variables. The
  script also loads the quantum variables, but are not used (so saving
  them is optional for this script, though you need to comment the
  relevant lines out).

- ``bilinears_ID.ipynb``: see the effect of the quantum modes to the
  matter content at the right-hand-side of the equations for K, KB
  (extrinsic curvature). You can tune the quantum numbers ``kmax``,
  ``lmax``, the Pauli-Villars mass ``mPV``, the convention for Planck
  mass, and more, and see how close the **initial** stress-energy
  tensor is to Minkowski, when the specific backreaction is
  considered. It might be useful to check this before running a
  semiclassical example with backreaction.

- ``calculate_UV_correlators.jl``:

- ``check_Hamiltonian_momentum_reduction_constraints.ipynb``:

- ``check_reduction_constraints_convergence.ipynb``:

- ``classical_bisection_echos_infalling_rmax_convergence.ipynb``:

- ``classical_bisection_echos_infalling_rmax_evol.ipynb``:

- ``classical_convergence.ipynb``:

- ``classical_norm_convergence.ipynb``:

- ``classical_norm_convergence.jl``:

- ``classical_norm_convergence_noisy.ipynb``:

- ``classical_pointwise_convergence.jl``:

- ``classical_scalar_evolution.ipynb``:

- ``correlators.ipynb``:

- ``evol_quantum_backreact_vs_no_backreact.ipynb``:

- ``evol_quantum.ipynb``:

- ``evol_quantum_standing_waves_mink_norms.jl``:

- ``evol_quantum_standing_waves_mink_plots.ipynb``:

- ``mink_backreact_consistency_check.ipynb``:

- ``quantum_setup_check_ID.ipynb``:

- ``UV_correlators.ipynb``:
