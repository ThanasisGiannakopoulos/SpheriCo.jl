.. _code_structure:

Code structure
===============

SpheriCo.jl is a Julia package which consists of three modules,
located under ``SpheriCo.jl/src/``. The main module is ``SpheriCo.jl``
and the submodules are ``classical.jl`` under
``SpheriCo.jl/src/classical/``, and ``quantum.jl`` under
``SpheriCo.jl/src/quantum/``. The former provides the infrastructure
to perform simulations in classical gravity, and the latter in
semiclassical. Below you can read more about the basic functionality
and structure of each module.

.. _code_structure_main_module:

Main module ``SpheriCo.jl``
---------------------------

This module creates and exports the functions and
`structures <https://docs.julialang.org/en/v1/base/base/#struct>`_
needed in both submodules to perform the simulations. Here is a list
of them, together with a short description:

- ``grid.jl``: provides the tools to build the spacial grid for the
  simulation. It involves the mutable structure ``Grid``, which reads
  the number of grid points ``Nr`` and the maximum value of the radius
  ``r_max``. This parameters may be given only at the beginning of the
  simulation, by the user, or altered during the evolution by the code
  automatically, if the user chooses to perform the simulations with
  an infalling ``r_max``. The spatial grid is built by the function
  ``System``.

- ``operators.jl``: creates the operators that approximate spatial
  derivatives (``Dr_FD2``, ``Drr_FD2``, ``Dr_SBP2``), provide
  numerical dissipation (``KO_FD2``, ``low_pass``), and perform
  projections of variables on a smaller grid
  (``project_on_smaller_grid``). The latter is needed for simulations
  with infalling ``r_max``.

- ``ghosts.jl``: provides the functions that populate the ghost point of
  the spatial grid (needed to calculate derivatives numerically at
  r=0).

- ``write_output.jl``: defines functions used to save different data
  during the simulation.

- ``utils.jl``: provides various functions needed during the simulations
  (apparent horizon finder ``find_AH``) or in post-processing (list
  all hdf5 files ``list_h5_files``, calculate the Ricci scalar
  ``Ricci_scalar``).

- ``time_integrators.jl``: functions that perform the timestep forward
  during the evolution. They assume the knowledge of the
  right-hand-side (denoted by ``F``) of the integrated
  equations. They are not exported, but explicitly included at the
  level of each of the submodules.

  .. _code_structure_submodule_classical:

Submodule ``classical.jl``
---------------------------

This submodule builds the function ``run_classical`` that performs the
simulation in classical gravity. It includes the following:

- ``parameters.jl``: creates a structure to hold all the parameters
  needed for the simulation. They are chosen by the user when an
  example is run. Not all parameters have default values and an error
  will appear if they are left empty.

- ``rhs.jl``: builds the right-hand-side of the evolved equations.

- ``../time_integrators.jl``: defined in the main module and included
  here by providing the relative path. They use the right-hand-side
  defined in the previous item.

- ``constraint.jl``: defines the function ``constraints`` that
  calculates the violation of the Hamiltonian and momentum
  constraints. It is used only in post-processing, but defined in two
  ways, such that it can also be used easier during the evolution.

- ``ID.jl``: constructs the initial data for the classical setup.

- ``utils.jl``: creates useful function for the simulation that
  provide output (``out``), or provide the option to exit the
  simulation when an apparent horizon forms (``AH_break``).

- ``evol.jl``: creates the function ``run_classical`` that performs
  the simulation for the classical setup. It needs the grid and
  parameters as inputs and uses the above functions, as well as those
  from the main module. The evolution mainly uses the `third order
  Adams-Bashforth method
  <https://en.wikipedia.org/wiki/Linear_multistep_method#Adams%E2%80%93Bashforth_methods>`_,
  apart from the first two steps, where the `fourth order Runge-Kutta
  method
  <https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#The_Runge%E2%80%93Kutta_method>`_
  is used.


  .. _code_structure_submodule_quantum:

Submodule ``quantum.jl``
---------------------------

This submodule builds the function ``run_quantum`` that performs the
simulation in the semiclassical setup. It includes the following:

- ``parameters.jl``: provides the structure with all the necessary
  parameters for the simulations. Some of them are the same as in the
  ``classical.jl`` submodule. Not all of the parameters have default
  values and if left empty an error will appear.

- ``../time_integrators.jl``: included from the main module (no
  exported functions).

- ``../classical/constraints.jl``: the classical Hamiltonian and
  momentum constraint violation, loaded from the classical submodule
  using the relative path.

- ``../classical/ID.jl``: the classical initial data, loaded from the
  classical submodule using the relative path.

- ``ID.jl``: the initial data for the quantum modes.

- ``bilinears.jl``: functions to calculate the bilinears, expressions
  necessary for the backreaction of quantum matter to geometry.

- ``correlators.jl``: functions to calculate correlators at equal
  simulation time, in space (2-point functions).

- ``rhs.jl``: the right-hand-side of the evolved equations. Includes
  both the classical equations (as in ``classical.jl``) and the
  equations of motion for the quantum modes (quantum matter). It has
  versions for reqularized and non-regularized quantum modes (state
  vectors of different dimensions), as well as for a semiclassical
  setup with and without backreaction (conditional statements).

- ``utils.jl``: various functions used during the evolution to save
  data, print output, or stop the simulation if an apparent horizon is
  formed.

- ``evol.jl``: creates the function ``run_quantum``, that performs the
  semiclassical simulation. It needs the grid and the model
  parameters. The time evolution is performed as in the classical
  case.
