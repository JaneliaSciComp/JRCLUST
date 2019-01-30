.. _tutorial:

Tutorial
========

.. _installing-jrclust:

Installing JRCLUST
------------------

JRCLUST is hosted on `GitHub <https://github.com/JaneliaSciComp/JRCLUST>`__.
If you'd like to test the latest development code, you can clone the repository to your computer.
If you want to stay on a release, head to the `releases page <https://github.com/JaneliaSciComp/JRCLUST/releases>`__ and download the latest release.

You may want to add something like the following to your `startup script <https://www.mathworks.com/help/matlab/ref/startup.html>`__:

.. code-block:: matlab

   addpath('/path/to/JRCLUST');

You may also need to recompile your CUDA codes if you're not on Windows.
Do this with

.. code-block:: matlab

   jrclust.CUDA.compileCUDA();

.. warning::

   You may get this error on recent versions of Ubuntu: ``Call to sgemv in CUBLAS failed with error status: CUBLAS_STATUS_EXECUTION_FAILED.``
   This is a `known issue <https://www.mathworks.com/matlabcentral/answers/437756-how-can-i-recompile-the-gpu-libraries>`__, but unfortunately the suggested workaround doesn't seem to work.

Setting up your config file
---------------------------

JRCLUST requires a configuration file specifying a number of :ref:`parameters <parameters>`.
How some of these parameters are used is explained in the :ref:`pipeline <pipeline>` section.
You may set this up with

.. code-block:: matlab

   jrc bootstrap

or

.. code-block:: matlab

   jrc bootstrap /path/to/metafile.meta

See the :ref:`bootstrap documentation <bootstrapping>` for details.

You should end up with a parameter file that looks something like this:

.. literalinclude:: ../../../default.prm
   :language: matlab

See the :ref:`parameters page <parameters` for details.
