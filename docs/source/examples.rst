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
   :code_block:`export OMP_NUM_THREADS=1` before running the script (see
   `here <https://github.com/JuliaLang/julia/issues/33409>`_ for a
   related discussion).

5. Run the script with :code_block:`julia run_classical.jl`. The data are
   saved in the file ``<root_dir>/data_<Nr>/``, where ``root_dir`` and
   ``Nr`` (number of points on the spatial grid) are parameters in
   ``run_classical.jl``.

.. _examples-simeclassical:

Semiclassical
------------
