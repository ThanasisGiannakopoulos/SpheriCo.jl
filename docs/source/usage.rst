.. _installation:

Installation
=============


Standard
------------

The package is not included in the Julia registry. After you download
Julia, you can download the package by typing ``julia`` in your
terminal, entering the package environment by pressing ``]`` and then
running

.. code-block:: console

   pkg> add https://github.com/ThanasisGiannakopoulos/SpheriCo.jl

Development
----------------

If you would like to make changes to SpheriCo.jl or contribute to its
development, you should first clone the repository. Then in your
terminal, go to your local directory of the repository and enter the
Julia REPL by typing ``julia``. Finally, enter the package environment
by pressing ``]`` and run


.. code-block:: console

   pkg> dev .

