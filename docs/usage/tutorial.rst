.. _tutorial:

Tutorial
========

.. _installing-jrclust:

Installing JRCLUST
------------------

JRCLUST is hosted on `GitHub <https://github.com/JaneliaSciComp/JRCLUST>`__.
If you'd like to test the latest development code, you can `clone the repository <https://help.github.com/articles/cloning-a-repository/>`__ to your computer.
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

.. literalinclude:: ../../default.prm
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

The preview GUI
---------------

A more advanced view on your raw traces is the preview GUI.
Invoke it like so:

.. code-block:: matlab

   jrc preview /path/to/your/configfile.prm

This will take a look at your recording and do some preliminary spike detection with the parameters you've specified.
(Spikes are circled in red below).
You can view the maximum site-to-site correlation and set a threshold for bad sites (i.e., :ref:`sites to ignore <ignoreSites>`) if they come below that threshold.
In the figure below, site 46 is poorly correlated with the other sites, so you might choose to ignore it.
You can also view the common average across sites at each time step, expressed in units of MAD, to set a :ref:`threshold <blankThresh>` for blanking out that period.
There are many other parameters to set from the preview GUI.
Take some time to explore these in the Edit menu, and see the effects they have by changing views in the View menu.
(You can see the :ref:`manual page <man-preview>` for details.)
Once you are satisfied with these parameters, you can select "Save to [your parameter file here]" from the File menu and start to detect spikes.

.. image:: /.static/plot-preview.PNG

.. _tut-detect:

Detecting spikes
----------------

Now it's time to detect spikes in your file.
This can be done with

.. code-block:: matlab

   jrc detect /path/to/your/configfile.prm

If you want to cluster your spikes immediately afterward, use

.. code-block:: matlab

   jrc detect-sort /path/to/your/configfile.prm

instead.

Depending on your choice of parameters, you should see something like the following output:

.. code-block:: matlab

   >> jrc detect test.prm
   Clearing GPU memory...done
   Processing load 1/1...
   Loading from file...done (0.06 s)
   Filtering spikes...  done (0.07) s
   Detecting spikes from each channel.
   ................................................................
   Detected 14883 spikes from 64 sites; took 0.3s.
   Merging spikes...  9891 spiking events found; took 0.1s
   Extracting features...done (2.39 s)
   File 1/1 took 8.3s (101.1 MB, 12.2 MB/s, x1.0 realtime)
   Detection completed in 8.30 seconds
   Saved spikesRaw to F:\Tests\JRCLUST\single\test_raw.jrc
   Saved spikesFilt to F:\Tests\JRCLUST\single\test_filt.jrc
   Saved spikeFeatures to F:\Tests\JRCLUST\single\test_features.jrc
   Saved results to F:\Tests\JRCLUST\single\test_res.mat

Your detection results will be exported to the workspace in the form of a struct called ``res`` for inspection.
See :ref:`io-files` for a description of the contents.
For a detailed description of the detect step, see :ref:`pipeline-detect`.

.. image:: /.static/detect.PNG

Clustering spikes
-----------------

Once you have detected your spikes, it's time to cluster them:

.. code-block:: matlab

   jrc sort /path/to/your/configfile.prm

You should see something like the following output:

.. code-block:: matlab

   >> jrc sort test.prm
   Using spikes detected on 31-Jan-2019 09:38:17
   Clearing GPU memory...done
   Computing rho
     ................................................................................................................................
     took 5.4s
   Computing delta
     ................................................................
     took 0.13s
   assigning clusters, nClusters:1177
   i: 1, n0 = 1869, i: 2, n0 = 1869, i: 3, n0 = 1869, i: 4, n0 = 1869, i: 5, n0 = 1869, i: 6, n0 = 1869, i: 7, n0 = 1869, i: 8, n0 = 1869, i: 9, n0 = 1869, i: 10, n0 = 1869,
     took 0.1s. Removed 1131 clusters having <30 spikes: ->46
   ..............................................Calculating cluster mean waveform.
     ............................................................................................
     took 0.1s
   Computing waveform correlation...  took 0.2s
   Computing self correlation
     ..............................................
     took 0.0s
   Calculating cluster mean waveform.
     ............................................................................................
     took 0.1s
   Computing waveform correlation...  took 0.1s
   Computing self correlation
     ..............................................
     took 0.0s
   Calculating cluster quality...
     took 0.0s
   Sorting completed in 6.34 seconds
   Saved spikesRaw to F:\Tests\JRCLUST\single\test_raw.jrc
   Saved spikesFilt to F:\Tests\JRCLUST\single\test_filt.jrc
   Saved spikeFeatures to F:\Tests\JRCLUST\single\test_features.jrc
   Saved results to F:\Tests\JRCLUST\single\test_res.mat

JRCLUST will tell you when your spikes were detected and perform :ref:`clustering <pipeline-sort>` and postprocessing.
As in the :ref:`detect step <tut-detect>`, JRCLUST will export the results structure to the workspace for inspection.

Curating your clustering
------------------------

You will now want to inspect the results of your clustering.
Do this with

.. code-block:: matlab

   jrc manual /path/to/your/configfile.prm

You will be greeted with the following screen:

.. image:: /.static/curate.PNG

Each of these figures contains a different view onto the data.
Here you will be able to annotate, delete, merge, or split clusters.
For a description of what each figure does and how to perform these operations, see the :ref:`pipeline-curate` section.
When you are satisfied with your clustering, you may save and exit.
