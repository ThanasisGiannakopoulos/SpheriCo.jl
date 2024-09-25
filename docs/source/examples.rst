.. _examples:

Examples
=============

Here you can find instructions on how to run the examples included
with the package, in the terminal. You could also run the examples in
Julia's REPL or in Jupyter notebooks, but instructions are not
provided and there may be differences. All examples are under
``SpheriCo.jl/examples/``.

.. _examples-classical:

Classical
------------

A single classical simulation can be run with
``run_classical.jl``. To do so follow the steps:

1. Go to the directory of the examples.

2. Open the script ``run_classical.jl`` with your favourite text editor.

3. Tune the parameters according to what you want to simulate. There
   is a description of each parameter within the script.

4. There is no parallelization option for the classical simulation, so
   running with 1 thread (the default option) is all that is
   possible. Sometimes, more threads might get activated during the
   simulations (not sure why, based on experience with simulations in
   Ubuntu 22.04). It might be useful to run the following command

   .. code-block:: console

      export OMP_NUM_THREADS=1

   before running the script (see
   `here <https://github.com/JuliaLang/julia/issues/33409>`_ for a
   related discussion).

5. Run the script with

   .. code-block:: console

      julia run_classical.jl

   The data are saved in the file ``<root_dir>/data_<Nr>/``, where
   ``root_dir`` and ``Nr`` (number of points on the spatial grid) are
   parameters in ``run_classical.jl``.

When you run the classical simulation with fixed ``rmax``, the output
should look like this

.. image:: ../images/run_classical.png
  :width: 800

The top row explains what is the output on each column:

1. Simulation speed (momentarily) in terms of steps/hour. It is
   initiated as zero and can appear as ``Inf`` (very fast, does not
   affect the simulation). This speed is smaller for the first two
   steps, since the time integrator there is the 4th order
   Runge-Kutta, which is slower than the 3rd order Adams-Bashforth
   used mainly.

2. Iteration (or step) of the simulation.

3. Time (in code units) of the simulation. The information from
   columns 1-3 can be used to make a rough estimate of how long a
   given simulation can take. This is more complicated with an
   infalling ``rmax`` setup, because in this case the timestep gets
   smaller as the simulation progresses.

4. The maximum value of the classical scalar field :math:`{\Phi}`.

5. The minimum value of :math:`{\Phi}`.

6. The maximum of the absolute value of the variable :math:`{\Theta}`.
   It is initiated to zero and is proportional to the Hamiltonian
   constraint violation.

7. The maximum of the absolute value of the variable :math:`{Z_r}`.
   It is also initiated to zero and is proportional to the momentum
   constraint violation.

8. The minimum of the lapse function :math:`{\alpha}`.

9. The maximum of :math:`{\alpha}`. The information on columns 7,8 can
   indicate whether and apparent horizon forms (if min is close to
   zero and max is close to 1).

   
If you run with infalling ``rmax``, the output should look like this

.. image:: ../images/run_classical_infalling_rmax.png
  :width: 800

In this case there are two extra columns (10,11) at the end:

10. The position of the apparent horizon ``r_AH``. If it is negative,
    there is no apparent horizon.

11. The position of the outer boundary ``rmax``.


.. _examples-simeclassical:

Semiclassical
------------
