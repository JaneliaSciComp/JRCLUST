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

You will be guided through a series of prompts in setting up your config file.

Selecting a meta file
~~~~~~~~~~~~~~~~~~~~~

If you specify a meta file with ``jrc bootstrap /path/to/metafile.meta``, then you will not be asked to select one;
you can skip to the next section.
After typing ``jrc bootstrap``, you should get a prompt asking if you have a meta file:

.. image:: /.static/bootstrap-meta-ask.PNG

If you select "Yes", you will be prompted to select one or more with a file dialog:

.. image:: /.static/bootstrap-meta-select.PNG

If you select "No" instead, you will be prompted to select one or more **raw recordings** with a similar dialog.

In either case, your working directory will be set to the directory containing the file or files you selected.

Selecting a probe file
~~~~~~~~~~~~~~~~~~~~~~

You will next be asked if you have a :ref:`probe file <io-probe>` you want to specify:

.. image:: /.static/bootstrap-probe-ask.PNG

If you select "Yes", you will be prompted to select exactly one probe file with a file dialog:

.. image:: /.static/bootstrap-probe-select.PNG

.. note::
   If you have a probe file in your working directory, then the dialog will search there first.
   If you do **not** have a probe file in your working directory, the dialog will search in the
   default location, JRCLUST/probes, as shown above.

Confirmation
~~~~~~~~~~~~

Next, you will be asked to confirm some salient data:

.. image:: /.static/bootstrap-confirm.PNG

Finishing up
~~~~~~~~~~~~

Once you are satisfied, JRCLUST will open up your parameter file for editing:

.. image:: /.static/bootstrap-finish.PNG

You should end up with a parameter file that looks something like this
(see the :ref:`parameters page <parameters>` for details):

.. literalinclude:: ../../../default.prm
   :language: matlab

Look over your parameters and ensure they are satisfactory.

Displaying your probe
---------------------

You might wish to plot your probe to ensure that it's modeled as expected.
This can be done with

.. code-block:: matlab

   jrc probe /path/to/your/configfile.prm

You will get a plot similar to this:

.. image:: /.static/plot-probe.PNG

In this example, we have a two-shank probe model, each shank having 32 active electrodes
for a total of 64 sites.
The shanks are spaced 250 μm apart, with the leftmost shank being considered the origin.
(Your probe will almost certainly look different.)

The labels on the patch plots are **site number**/**channel number**.
This is handy for visualizing the channel-site correspondence.
You can use this graphic to verify that the sites you gave to JRCLUST are in their proper positions.
You can zoom in and out with the scroll wheel to increase or decrease resolution as necessary.

Close this plot and proceed to the next section.

Showing raw traces
------------------

You will probably want to take a look at your raw traces.
Do this with

.. code-block:: matlab

   jrc traces /path/to/your/configfile.prm

You will get an interactive plot that looks like this:

.. image:: /.static/plot-traces.PNG

Here you can page through your raw data, adjusting your scale and applying various filters.
Later, after you have detected and clustered your spikes, you can highlight spikes like so:

.. image:: /.static/plot-traces-detections.PNG

You can also plot the `power spectral density <https://en.wikipedia.org/wiki/Spectral_density>`__ by hitting the ``P`` key in this plot.
You may select a single site or show the mean power vs. frequency as shown below:

.. image:: /.static/plot-psd.PNG
