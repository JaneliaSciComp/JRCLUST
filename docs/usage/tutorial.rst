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
You may select from among SpikeGLX .bin/.dat files, or Intan .rhd traditional-format files.

.. image:: /.static/bootstrap-intan.PNG

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
  2019-08-30 14:04:09 Clearing GPU memory...
  2019-08-30 14:04:09 GPU memory cleared (took 0.16 s)
  2019-08-30 14:04:09 Detecting spikes in recordings...
  2019-08-30 14:04:09 Processing file F:\Misc\Tests\JRCLUST\single\test.bin (1/1)...
  2019-08-30 14:04:10 Processing load 1/1...
  2019-08-30 14:04:10 Filtering samples...
  2019-08-30 14:04:10 Finished filtering samples (took 0.10 s)
  2019-08-30 14:04:11 Detecting spikes from each site...
  2019-08-30 14:04:11 Detected 19 spikes on site 1
  2019-08-30 14:04:11 Detected 17 spikes on site 2
  2019-08-30 14:04:11 Detected 69 spikes on site 3
  2019-08-30 14:04:11 Detected 127 spikes on site 4
  2019-08-30 14:04:11 Detected 29 spikes on site 5
  2019-08-30 14:04:11 Detected 22 spikes on site 6
  2019-08-30 14:04:11 Detected 51 spikes on site 7
  2019-08-30 14:04:11 Detected 61 spikes on site 8
  2019-08-30 14:04:11 Detected 52 spikes on site 9
  2019-08-30 14:04:11 Detected 31 spikes on site 10
  2019-08-30 14:04:11 Detected 113 spikes on site 11
  2019-08-30 14:04:11 Detected 139 spikes on site 12
  2019-08-30 14:04:11 Detected 61 spikes on site 13
  2019-08-30 14:04:11 Detected 122 spikes on site 14
  2019-08-30 14:04:11 Detected 274 spikes on site 15
  2019-08-30 14:04:11 Detected 1156 spikes on site 16
  2019-08-30 14:04:11 Detected 223 spikes on site 17
  2019-08-30 14:04:11 Detected 83 spikes on site 18
  2019-08-30 14:04:11 Detected 167 spikes on site 19
  2019-08-30 14:04:11 Detected 362 spikes on site 20
  2019-08-30 14:04:11 Detected 528 spikes on site 21
  2019-08-30 14:04:11 Detected 248 spikes on site 22
  2019-08-30 14:04:11 Detected 282 spikes on site 23
  2019-08-30 14:04:11 Detected 331 spikes on site 24
  2019-08-30 14:04:11 Detected 213 spikes on site 25
  2019-08-30 14:04:11 Detected 207 spikes on site 26
  2019-08-30 14:04:11 Detected 112 spikes on site 27
  2019-08-30 14:04:11 Detected 426 spikes on site 28
  2019-08-30 14:04:11 Detected 368 spikes on site 29
  2019-08-30 14:04:11 Detected 227 spikes on site 30
  2019-08-30 14:04:11 Detected 323 spikes on site 31
  2019-08-30 14:04:11 Detected 404 spikes on site 32
  2019-08-30 14:04:11 Detected 2047 spikes on site 33
  2019-08-30 14:04:11 Detected 150 spikes on site 34
  2019-08-30 14:04:11 Detected 69 spikes on site 35
  2019-08-30 14:04:11 Detected 68 spikes on site 36
  2019-08-30 14:04:11 Detected 121 spikes on site 37
  2019-08-30 14:04:11 Detected 97 spikes on site 38
  2019-08-30 14:04:11 Detected 102 spikes on site 39
  2019-08-30 14:04:11 Detected 89 spikes on site 40
  2019-08-30 14:04:11 Detected 114 spikes on site 41
  2019-08-30 14:04:11 Detected 215 spikes on site 42
  2019-08-30 14:04:11 Detected 251 spikes on site 43
  2019-08-30 14:04:11 Detected 457 spikes on site 44
  2019-08-30 14:04:11 Detected 645 spikes on site 45
  2019-08-30 14:04:11 Detected 182 spikes on site 46
  2019-08-30 14:04:11 Detected 278 spikes on site 47
  2019-08-30 14:04:11 Detected 408 spikes on site 48
  2019-08-30 14:04:11 Detected 337 spikes on site 49
  2019-08-30 14:04:11 Detected 347 spikes on site 50
  2019-08-30 14:04:11 Detected 168 spikes on site 51
  2019-08-30 14:04:11 Detected 118 spikes on site 52
  2019-08-30 14:04:11 Detected 196 spikes on site 53
  2019-08-30 14:04:11 Detected 185 spikes on site 54
  2019-08-30 14:04:11 Detected 105 spikes on site 55
  2019-08-30 14:04:11 Detected 275 spikes on site 56
  2019-08-30 14:04:11 Detected 121 spikes on site 57
  2019-08-30 14:04:11 Detected 116 spikes on site 58
  2019-08-30 14:04:11 Detected 250 spikes on site 59
  2019-08-30 14:04:11 Detected 232 spikes on site 60
  2019-08-30 14:04:11 Detected 62 spikes on site 61
  2019-08-30 14:04:11 Detected 139 spikes on site 62
  2019-08-30 14:04:11 Detected 51 spikes on site 63
  2019-08-30 14:04:11 Detected 41 spikes on site 64
  2019-08-30 14:04:11 Finished detecting spikes (took 0.15 s)
  2019-08-30 14:04:11 Merging duplicate spiking events...
  2019-08-30 14:04:11 9891 spiking events found (took 0.10 s)
  2019-08-30 14:04:11 Finished load 1/1 (took 1.46 s)
  2019-08-30 14:04:11 Finished processing file F:\Misc\Tests\JRCLUST\single\test.bin (1/1) (took 2.04 s)
  2019-08-30 14:04:11 Extracting features...
  2019-08-30 14:04:12 Finished extracting features (took 0.31 s)
  2019-08-30 14:04:12 Finished detecting (took 2.50 s)
  2019-08-30 14:04:15 Saving results to F:\Misc\Tests\JRCLUST\single\test_res.mat...
  2019-08-30 14:04:15 Results saved to F:\Misc\Tests\JRCLUST\single\test_res.mat (took 0.02 s)

  ====DETECTION SUMMARY====
  Detection completed in 2.50 s
  Spike count: 9891
  Spike counts per site: min 8 (site 2), max 2036 (site 33), median 85

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
  2019-08-30 14:06:37 Loading F:\Misc\Tests\JRCLUST\single\test_res.mat...
  2019-08-30 14:06:37 Finished loading F:\Misc\Tests\JRCLUST\single\test_res.mat (took 0.01 s)
  2019-08-30 14:06:37 Loading F:\Misc\Tests\JRCLUST\single\test_raw.jrc...
  2019-08-30 14:06:37 Finished loading F:\Misc\Tests\JRCLUST\single\test_raw.jrc (took 0.01 s)
  2019-08-30 14:06:37 Loading F:\Misc\Tests\JRCLUST\single\test_filt.jrc...
  2019-08-30 14:06:37 Finished loading F:\Misc\Tests\JRCLUST\single\test_filt.jrc (took 0.01 s)
  2019-08-30 14:06:37 Loading F:\Misc\Tests\JRCLUST\single\test_features.jrc...
  2019-08-30 14:06:37 Finished loading F:\Misc\Tests\JRCLUST\single\test_features.jrc (took 0.00 s)
  2019-08-30 14:06:37 Using spikes detected on 30-Aug-2019 14:04:12
  2019-08-30 14:06:37 Clearing GPU memory...
  2019-08-30 14:06:37 GPU memory cleared (took 0.16 s)
  2019-08-30 14:06:37 Sorting detected spikes...
  2019-08-30 14:06:37 Computing rho...
  2019-08-30 14:06:38 Site 1: rho cutoff, 731869.75; average rho: 0.00016 (15 spikes)
  2019-08-30 14:06:38 Site 2: rho cutoff, 1367942.13; average rho: 0.00036 (8 spikes)
  2019-08-30 14:06:38 Site 3: rho cutoff, 371085.19; average rho: 0.00059 (44 spikes)
  2019-08-30 14:06:38 Site 4: rho cutoff, 454989.97; average rho: 0.00098 (110 spikes)
  2019-08-30 14:06:38 Site 5: rho cutoff, 818252.50; average rho: 0.00114 (21 spikes)
  2019-08-30 14:06:38 Site 6: rho cutoff, 1096531.38; average rho: 0.00129 (8 spikes)
  2019-08-30 14:06:38 Site 7: rho cutoff, 575589.00; average rho: 0.00151 (34 spikes)
  2019-08-30 14:06:38 Site 8: rho cutoff, 481210.56; average rho: 0.00175 (29 spikes)
  2019-08-30 14:06:38 Site 9: rho cutoff, 636892.19; average rho: 0.00193 (40 spikes)
  2019-08-30 14:06:38 Site 10: rho cutoff, 2221163.75; average rho: 0.00211 (12 spikes)
  2019-08-30 14:06:38 Site 11: rho cutoff, 692781.94; average rho: 0.00248 (91 spikes)
  2019-08-30 14:06:38 Site 12: rho cutoff, 503940.38; average rho: 0.00293 (103 spikes)
  2019-08-30 14:06:38 Site 13: rho cutoff, 745752.50; average rho: 0.00311 (39 spikes)
  2019-08-30 14:06:38 Site 14: rho cutoff, 1094450.50; average rho: 0.00332 (43 spikes)
  2019-08-30 14:06:38 Site 15: rho cutoff, 855698.94; average rho: 0.00402 (222 spikes)
  2019-08-30 14:06:38 Site 16: rho cutoff, 308155.81; average rho: 0.00683 (1007 spikes)
  2019-08-30 14:06:38 Site 17: rho cutoff, 829652.13; average rho: 0.00736 (161 spikes)
  2019-08-30 14:06:38 Site 18: rho cutoff, 783599.13; average rho: 0.00758 (32 spikes)
  2019-08-30 14:06:38 Site 19: rho cutoff, 698739.75; average rho: 0.00802 (85 spikes)
  2019-08-30 14:06:38 Site 20: rho cutoff, 778769.88; average rho: 0.00854 (168 spikes)
  2019-08-30 14:06:38 Site 21: rho cutoff, 894751.50; average rho: 0.00948 (369 spikes)
  2019-08-30 14:06:38 Site 22: rho cutoff, 895366.63; average rho: 0.00986 (113 spikes)
  2019-08-30 14:06:38 Site 23: rho cutoff, 811321.75; average rho: 0.01045 (141 spikes)
  2019-08-30 14:06:38 Site 24: rho cutoff, 961806.00; average rho: 0.01101 (211 spikes)
  2019-08-30 14:06:38 Site 25: rho cutoff, 917075.13; average rho: 0.01133 (74 spikes)
  2019-08-30 14:06:38 Site 26: rho cutoff, 685997.25; average rho: 0.01184 (150 spikes)
  2019-08-30 14:06:38 Site 27: rho cutoff, 1064097.88; average rho: 0.01209 (47 spikes)
  2019-08-30 14:06:38 Site 28: rho cutoff, 539079.56; average rho: 0.01316 (250 spikes)
  2019-08-30 14:06:38 Site 29: rho cutoff, 1070611.50; average rho: 0.01401 (261 spikes)
  2019-08-30 14:06:38 Site 30: rho cutoff, 827772.50; average rho: 0.01425 (59 spikes)
  2019-08-30 14:06:38 Site 31: rho cutoff, 1146530.88; average rho: 0.01492 (153 spikes)
  2019-08-30 14:06:38 Site 32: rho cutoff, 970855.06; average rho: 0.01557 (241 spikes)
  2019-08-30 14:06:38 Site 33: rho cutoff, 251356.16; average rho: 0.02157 (2036 spikes)
  2019-08-30 14:06:38 Site 34: rho cutoff, 260651.09; average rho: 0.02213 (139 spikes)
  2019-08-30 14:06:38 Site 35: rho cutoff, 714987.13; average rho: 0.02235 (45 spikes)
  2019-08-30 14:06:38 Site 36: rho cutoff, 631074.38; average rho: 0.02262 (52 spikes)
  2019-08-30 14:06:38 Site 37: rho cutoff, 453064.25; average rho: 0.02289 (60 spikes)
  2019-08-30 14:06:38 Site 38: rho cutoff, 805470.44; average rho: 0.02323 (71 spikes)
  2019-08-30 14:06:38 Site 39: rho cutoff, 660921.75; average rho: 0.02345 (45 spikes)
  2019-08-30 14:06:38 Site 40: rho cutoff, 685622.75; average rho: 0.02372 (39 spikes)
  2019-08-30 14:06:38 Site 41: rho cutoff, 1017240.00; average rho: 0.02414 (79 spikes)
  2019-08-30 14:06:38 Site 42: rho cutoff, 724510.81; average rho: 0.02464 (129 spikes)
  2019-08-30 14:06:39 Site 43: rho cutoff, 1318095.13; average rho: 0.02502 (85 spikes)
  2019-08-30 14:06:39 Site 44: rho cutoff, 810584.25; average rho: 0.02592 (312 spikes)
  2019-08-30 14:06:39 Site 45: rho cutoff, 913448.50; average rho: 0.02711 (374 spikes)
  2019-08-30 14:06:39 Site 46: rho cutoff, 928920.94; average rho: 0.02750 (67 spikes)
  2019-08-30 14:06:39 Site 47: rho cutoff, 1020212.25; average rho: 0.02830 (236 spikes)
  2019-08-30 14:06:39 Site 48: rho cutoff, 465102.56; average rho: 0.02917 (257 spikes)
  2019-08-30 14:06:39 Site 49: rho cutoff, 804615.63; average rho: 0.02972 (188 spikes)
  2019-08-30 14:06:39 Site 50: rho cutoff, 559785.13; average rho: 0.03033 (187 spikes)
  2019-08-30 14:06:39 Site 51: rho cutoff, 557701.31; average rho: 0.03050 (16 spikes)
  2019-08-30 14:06:39 Site 52: rho cutoff, 371283.31; average rho: 0.03069 (39 spikes)
  2019-08-30 14:06:39 Site 53: rho cutoff, 542745.94; average rho: 0.03103 (94 spikes)
  2019-08-30 14:06:39 Site 54: rho cutoff, 498734.75; average rho: 0.03148 (137 spikes)
  2019-08-30 14:06:39 Site 55: rho cutoff, 1296482.13; average rho: 0.03166 (20 spikes)
  2019-08-30 14:06:39 Site 56: rho cutoff, 377575.56; average rho: 0.03241 (241 spikes)
  2019-08-30 14:06:39 Site 57: rho cutoff, 347963.09; average rho: 0.03268 (59 spikes)
  2019-08-30 14:06:39 Site 58: rho cutoff, 849382.00; average rho: 0.03297 (58 spikes)
  2019-08-30 14:06:39 Site 59: rho cutoff, 531617.50; average rho: 0.03337 (102 spikes)
  2019-08-30 14:06:39 Site 60: rho cutoff, 354038.84; average rho: 0.03401 (177 spikes)
  2019-08-30 14:06:39 Site 61: rho cutoff, 357702.53; average rho: 0.03422 (26 spikes)
  2019-08-30 14:06:39 Site 62: rho cutoff, 881595.94; average rho: 0.03482 (123 spikes)
  2019-08-30 14:06:39 Site 63: rho cutoff, 457689.47; average rho: 0.03507 (36 spikes)
  2019-08-30 14:06:39 Site 64: rho cutoff, 699817.06; average rho: 0.03532 (21 spikes)
  2019-08-30 14:06:39 Finished computing rho (took 1.27 s)
  2019-08-30 14:06:39 Computing delta...
  2019-08-30 14:06:39 Site 1: median delta: 1.27137 (15 spikes)
  2019-08-30 14:06:39 Site 2: median delta: 1.66295 (8 spikes)
  2019-08-30 14:06:39 Site 3: median delta: 1.00265 (44 spikes)
  2019-08-30 14:06:39 Site 4: median delta: 0.85708 (110 spikes)
  2019-08-30 14:06:39 Site 5: median delta: 1.58162 (21 spikes)
  2019-08-30 14:06:39 Site 6: median delta: 1.52481 (8 spikes)
  2019-08-30 14:06:39 Site 7: median delta: 1.24786 (34 spikes)
  2019-08-30 14:06:39 Site 8: median delta: 1.22150 (29 spikes)
  2019-08-30 14:06:39 Site 9: median delta: 1.08609 (40 spikes)
  2019-08-30 14:06:39 Site 10: median delta: 1.10892 (12 spikes)
  2019-08-30 14:06:39 Site 11: median delta: 0.87843 (91 spikes)
  2019-08-30 14:06:39 Site 12: median delta: 0.86909 (103 spikes)
  2019-08-30 14:06:39 Site 13: median delta: 1.10214 (39 spikes)
  2019-08-30 14:06:39 Site 14: median delta: 1.15297 (43 spikes)
  2019-08-30 14:06:39 Site 15: median delta: 0.73530 (222 spikes)
  2019-08-30 14:06:39 Site 16: median delta: 0.55511 (1007 spikes)
  2019-08-30 14:06:39 Site 17: median delta: 0.75445 (161 spikes)
  2019-08-30 14:06:39 Site 18: median delta: 1.07030 (32 spikes)
  2019-08-30 14:06:39 Site 19: median delta: 0.90345 (85 spikes)
  2019-08-30 14:06:39 Site 20: median delta: 0.80473 (168 spikes)
  2019-08-30 14:06:39 Site 21: median delta: 0.66820 (369 spikes)
  2019-08-30 14:06:39 Site 22: median delta: 0.91149 (113 spikes)
  2019-08-30 14:06:39 Site 23: median delta: 0.76094 (141 spikes)
  2019-08-30 14:06:39 Site 24: median delta: 0.75902 (211 spikes)
  2019-08-30 14:06:39 Site 25: median delta: 0.97118 (74 spikes)
  2019-08-30 14:06:39 Site 26: median delta: 0.81842 (150 spikes)
  2019-08-30 14:06:39 Site 27: median delta: 0.95115 (47 spikes)
  2019-08-30 14:06:39 Site 28: median delta: 0.72237 (250 spikes)
  2019-08-30 14:06:39 Site 29: median delta: 0.69283 (261 spikes)
  2019-08-30 14:06:39 Site 30: median delta: 1.03695 (59 spikes)
  2019-08-30 14:06:39 Site 31: median delta: 0.75059 (153 spikes)
  2019-08-30 14:06:39 Site 32: median delta: 0.74415 (241 spikes)
  2019-08-30 14:06:39 Site 33: median delta: 0.47524 (2036 spikes)
  2019-08-30 14:06:39 Site 34: median delta: 0.83504 (139 spikes)
  2019-08-30 14:06:39 Site 35: median delta: 1.03859 (45 spikes)
  2019-08-30 14:06:39 Site 36: median delta: 0.97808 (52 spikes)
  2019-08-30 14:06:39 Site 37: median delta: 0.98207 (60 spikes)
  2019-08-30 14:06:39 Site 38: median delta: 0.94462 (71 spikes)
  2019-08-30 14:06:39 Site 39: median delta: 1.00673 (45 spikes)
  2019-08-30 14:06:39 Site 40: median delta: 0.99713 (39 spikes)
  2019-08-30 14:06:39 Site 41: median delta: 0.87949 (79 spikes)
  2019-08-30 14:06:39 Site 42: median delta: 0.84302 (129 spikes)
  2019-08-30 14:06:39 Site 43: median delta: 0.90183 (85 spikes)
  2019-08-30 14:06:39 Site 44: median delta: 0.68221 (312 spikes)
  2019-08-30 14:06:39 Site 45: median delta: 0.66591 (374 spikes)
  2019-08-30 14:06:39 Site 46: median delta: 0.92985 (67 spikes)
  2019-08-30 14:06:39 Site 47: median delta: 0.70058 (236 spikes)
  2019-08-30 14:06:39 Site 48: median delta: 0.73736 (257 spikes)
  2019-08-30 14:06:39 Site 49: median delta: 0.76567 (188 spikes)
  2019-08-30 14:06:39 Site 50: median delta: 0.80459 (187 spikes)
  2019-08-30 14:06:39 Site 51: median delta: 1.46500 (16 spikes)
  2019-08-30 14:06:39 Site 52: median delta: 1.09898 (39 spikes)
  2019-08-30 14:06:39 Site 53: median delta: 0.92572 (94 spikes)
  2019-08-30 14:06:39 Site 54: median delta: 0.81506 (137 spikes)
  2019-08-30 14:06:39 Site 55: median delta: 1.22181 (20 spikes)
  2019-08-30 14:06:39 Site 56: median delta: 0.73346 (241 spikes)
  2019-08-30 14:06:39 Site 57: median delta: 1.03142 (59 spikes)
  2019-08-30 14:06:39 Site 58: median delta: 0.98693 (58 spikes)
  2019-08-30 14:06:39 Site 59: median delta: 0.86694 (102 spikes)
  2019-08-30 14:06:39 Site 60: median delta: 0.83345 (177 spikes)
  2019-08-30 14:06:39 Site 61: median delta: 1.02984 (26 spikes)
  2019-08-30 14:06:39 Site 62: median delta: 0.77846 (123 spikes)
  2019-08-30 14:06:39 Site 63: median delta: 1.04947 (36 spikes)
  2019-08-30 14:06:39 Site 64: median delta: 1.07382 (21 spikes)
  2019-08-30 14:06:39 Finished computing delta (took 0.27 s)
  2019-08-30 14:06:39 Assigning clusters (nClusters: 528)...
  2019-08-30 14:06:39 9363/9891 spikes unassigned, 709 can be assigned
  2019-08-30 14:06:39 8654/9891 spikes unassigned, 1015 can be assigned
  2019-08-30 14:06:39 7639/9891 spikes unassigned, 1144 can be assigned
  2019-08-30 14:06:39 6495/9891 spikes unassigned, 1065 can be assigned
  2019-08-30 14:06:39 5430/9891 spikes unassigned, 977 can be assigned
  2019-08-30 14:06:39 4453/9891 spikes unassigned, 866 can be assigned
  2019-08-30 14:06:39 3587/9891 spikes unassigned, 740 can be assigned
  2019-08-30 14:06:39 2847/9891 spikes unassigned, 620 can be assigned
  2019-08-30 14:06:39 2227/9891 spikes unassigned, 519 can be assigned
  2019-08-30 14:06:39 1708/9891 spikes unassigned, 438 can be assigned
  2019-08-30 14:06:39 1270/9891 spikes unassigned, 333 can be assigned
  2019-08-30 14:06:39 937/9891 spikes unassigned, 247 can be assigned
  2019-08-30 14:06:39 690/9891 spikes unassigned, 218 can be assigned
  2019-08-30 14:06:39 472/9891 spikes unassigned, 148 can be assigned
  2019-08-30 14:06:39 324/9891 spikes unassigned, 103 can be assigned
  2019-08-30 14:06:39 221/9891 spikes unassigned, 58 can be assigned
  2019-08-30 14:06:39 163/9891 spikes unassigned, 34 can be assigned
  2019-08-30 14:06:39 129/9891 spikes unassigned, 27 can be assigned
  2019-08-30 14:06:39 102/9891 spikes unassigned, 28 can be assigned
  2019-08-30 14:06:39 74/9891 spikes unassigned, 21 can be assigned
  2019-08-30 14:06:39 53/9891 spikes unassigned, 15 can be assigned
  2019-08-30 14:06:39 38/9891 spikes unassigned, 18 can be assigned
  2019-08-30 14:06:39 20/9891 spikes unassigned, 12 can be assigned
  2019-08-30 14:06:39 8/9891 spikes unassigned, 4 can be assigned
  2019-08-30 14:06:39 4/9891 spikes unassigned, 1 can be assigned
  2019-08-30 14:06:39 3/9891 spikes unassigned, 2 can be assigned
  2019-08-30 14:06:39 1/9891 spikes unassigned, 1 can be assigned
  2019-08-30 14:06:39 9826/9891 spikes unassigned, 272 can be assigned
  2019-08-30 14:06:39 9554/9891 spikes unassigned, 589 can be assigned
  2019-08-30 14:06:39 8965/9891 spikes unassigned, 867 can be assigned
  2019-08-30 14:06:39 8098/9891 spikes unassigned, 958 can be assigned
  2019-08-30 14:06:39 7140/9891 spikes unassigned, 1008 can be assigned
  2019-08-30 14:06:39 6132/9891 spikes unassigned, 976 can be assigned
  2019-08-30 14:06:39 5156/9891 spikes unassigned, 872 can be assigned
  2019-08-30 14:06:39 4284/9891 spikes unassigned, 726 can be assigned
  2019-08-30 14:06:39 3558/9891 spikes unassigned, 618 can be assigned
  2019-08-30 14:06:39 2940/9891 spikes unassigned, 529 can be assigned
  2019-08-30 14:06:39 2411/9891 spikes unassigned, 416 can be assigned
  2019-08-30 14:06:39 1995/9891 spikes unassigned, 325 can be assigned
  2019-08-30 14:06:39 1670/9891 spikes unassigned, 284 can be assigned
  2019-08-30 14:06:39 1386/9891 spikes unassigned, 185 can be assigned
  2019-08-30 14:06:39 1201/9891 spikes unassigned, 122 can be assigned
  2019-08-30 14:06:39 1079/9891 spikes unassigned, 69 can be assigned
  2019-08-30 14:06:39 1010/9891 spikes unassigned, 39 can be assigned
  2019-08-30 14:06:39 971/9891 spikes unassigned, 27 can be assigned
  2019-08-30 14:06:39 944/9891 spikes unassigned, 28 can be assigned
  2019-08-30 14:06:39 916/9891 spikes unassigned, 22 can be assigned
  2019-08-30 14:06:39 894/9891 spikes unassigned, 17 can be assigned
  2019-08-30 14:06:39 877/9891 spikes unassigned, 19 can be assigned
  2019-08-30 14:06:39 858/9891 spikes unassigned, 13 can be assigned
  2019-08-30 14:06:39 845/9891 spikes unassigned, 4 can be assigned
  2019-08-30 14:06:39 841/9891 spikes unassigned, 1 can be assigned
  2019-08-30 14:06:39 840/9891 spikes unassigned, 3 can be assigned
  2019-08-30 14:06:39 837/9891 spikes unassigned, 1 can be assigned
  2019-08-30 14:06:39 836/9891 spikes unassigned, 1 can be assigned
  2019-08-30 14:06:39 Finished initial assignment (65 clusters) (took 0.11 s)
  2019-08-30 14:06:39 Computing cluster mean waveforms...
  2019-08-30 14:06:39 Finished computing cluster mean waveforms (took 0.17 s)
  2019-08-30 14:06:39 Computing waveform correlation...
  2019-08-30 14:06:40 Finished computing waveform correlation (took 0.19 s)
  2019-08-30 14:06:40 Computing cluster self-similarity...
  2019-08-30 14:06:40 Finished computing cluster self-similarity (took 0.05 s)
  2019-08-30 14:06:40 Computing cluster quality scores...
  2019-08-30 14:06:40 Finished computing cluster quality scores (took 0.06 s)
  2019-08-30 14:06:40 Merging clusters by similarity...
  2019-08-30 14:06:40 No clusters to merge (took 0.00 s)
  2019-08-30 14:06:40 Finished sorting (took 2.52 s)
  2019-08-30 14:06:42 Saving results to F:\Misc\Tests\JRCLUST\single\test_res.mat...
  2019-08-30 14:06:43 Results saved to F:\Misc\Tests\JRCLUST\single\test_res.mat (took 0.13 s)

  ====SORTING SUMMARY====
  Sorting completed in 2.51 s
  Clusters: 65 (no merges)
  Spike count per cluster: min 31 (cluster 24), max 2036 (cluster 37), median 80
  Site count per cluster: min 1 (cluster 1), max 1 (cluster 1), median 1

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
