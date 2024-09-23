.. _code_structure:

Code structure
===============

SpheriCo.jl is a Julia package which consists of three modules,
located under ``SpheriCo.jl/src/``. The main module is
``SpheriCo.jl`` and the submodules are ``classical.jl`` under
``SpheriCo.jl/src/classical/``, and ``quantum.jl`` under
``SpheriCo.jl/src/quantum/``. The former submodule provides the
infrastructure to perform simulations in classical gravity, and the
latter in semiclassical. Below you can read more about the basic
functionality and structure of each module.

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
  simulation, by the user, or updated (changed) during the evolution
  by the code automatically, if the user chooses to perform the
  simulations with an infalling ``r_max``. The spatial grid is built
  by the function ``System``.

- ``operators.jl``: creates and exports the operators that approximate
  spatial derivatives (``Dr_FD2``, ``Drr_FD2``, ``Dr_SBP2``), provide
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
  right-hand-side (denoted by ``F`` + more test) of the integrated
  equations. They are not exported, but explicitly included at the
  level of each of the submodules. To be exported at the level of the
  main module, their structure needs to be modified appropriately.
