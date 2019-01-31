.. _io-files:

Input and output files
======================

Replace "mysession", "myprobe", and "myrecording" as appropriate.

Input files
-----------

+----------------------+---------------------+------------------+
| Filename             | Content             | Format           |
+======================+=====================+==================+
| mysession.prm        | :ref:`io-config`    | Plain text       |
|                      |                     | (MATLAB          |
|                      |                     | ``eval``-uable)  |
+----------------------+---------------------+------------------+
| myprobe.prb          | :ref:`io-probe`     | Plain text       |
|                      |                     | (MATLAB          |
|                      |                     | ``eval``-uable)  |
+----------------------+---------------------+------------------+
| myrecording.bin      | :ref:`io-raw`       | Binary (16-bit   |
|                      |                     | signed integer)  |
+----------------------+---------------------+------------------+
| myrecording.meta     | Meta file for       | Binary file      |
|                      | your recording      |                  |
|                      | (see :ref:`io-raw`) |                  |
+----------------------+---------------------+------------------+
| mysession\_trial.mat | :ref:`io-trial`     | MATLAB data or   |
|                      |                     | CSV              |
+----------------------+---------------------+------------------+

.. _io-config:

Config file
~~~~~~~~~~~

A MATLAB ``eval``-uable file containing your parameters.
See, e.g., :ref:`bootstrapping` and :ref:`parameters`.

.. _io-probe:

Probe file
~~~~~~~~~~

A MATLAB ``eval``-uable file containing parameters for your probe.
(Probe files for various configurations are stored in JRCLUST/probes/.)

.. note::

   Probe files are now **optional**.
   Probe file support is for backwards compatibility.
   Going forward, you can specify your probe parameters in your config file directly.

Your probe may be visualized with the ``probe`` :ref:`command <usage-cli>`.

Example probe file
^^^^^^^^^^^^^^^^^^

.. code-block:: matlab

    % Order of the probe sites in the recording file
    channels = 1 + [40  103 39  104 41  102 38  105 42  101 37  106 43  100 36  107 44  99  35  108 45  98  34  109 46  97  33  110 47  96  32  111 48  127 63  112 49  126 62  113 50  125 61  114 51  124 60  115 52  116 59  123 53  117 58  122 54  118 57  121 55  119 56  120 8   71  7   72  9   70  6   73  10  74  5   69  11  75  4   68 12   76  3   67  13  77  2   66  14  78  1   65  15  79  0   64 16   80  31  95  17  81  30  94  18  82  29  93  19  83  28  92     20   84  27  91  21  85  26  90  22  86  25  89  23  87  24  88];

    % Site location in micrometers (x and y)
    geometry = zeros(128, 2);    % Create a 128 x 2 matrix filled with zero
    geometry(1:2:end,1) = 0;     % Set the x-coordinates of the left-edge sites
    geometry(2:2:end,1) = 28;    % Set the x-coordinates of the right-edge sites
    geometry(1:2:end,2) = 20*(0:63);    % Set the y-coordinates of the left-edge sites
    geometry(2:2:end,2) = geometry(1:2:end,2); % Set the y-coordinates of the right edge sites by copying the left-edge site values

    % Reference sites are being excluded
    ref_sites = [1 18 33 50 65 82 97 114]; % Specify the site numbers to exclude
    channels(ref_sites) = [];  % Delete reference sites from the channels vector
    geometry(ref_sites,:) = []; % Delete reference sites from the geometry matrix

    % Recording contact pad size in micrometers. Height x width
    pad = [12 12];

    % Default prm
    maxSite = 4.5;    % Used to calculate the number of sites to group spikes (nSites_spk = 1 + maxSite*2)
    um_per_pix = 20;  % Vertical site center-to-center spacing (used for display)

    % Shank information
    shank = ones(size(channels));  % All sites belong to shank 1 (single shank recording)

Tags (variable names)
^^^^^^^^^^^^^^^^^^^^^

- ``channels``: (1 x ``nSites``), where ``nSites`` is the number of sites to load.
  The :math:`i` th entry in ``channels`` denotes the channel (the row in the :math:`n_{\text{channels}} \times n_{\text{samples}}` matrix) in the raw recording containing
  data from site :math:`i`.
  In the example above, channel 41 in the raw recording corresponds to site 1.
  This corresponds to the :ref:`siteMap` parameter.
- ``geometry``: (``nSites`` x 2).
  The location of each site, in μm.
  The first column corresponds to the horizontal or ``x`` dimension,
  and the second column corresponds to the vertical or ``y`` dimension (parallel to the probe shank).
  This corresponds to the :ref:`siteLoc` parameter.
- ``pad``: (1 x 2).
  Dimensions (height, width) of the rectangular recording site, in μm.
  This corresponds to the :ref:`probePad` parameter.
- ``ref_sites`` (optional): (1 x ``k``), where ``k`` is the number of reference sites.
  Indices of reference or disconnected sites to ignore.
  Alternatively, you can leave these sites out of your :ref:`siteMap` entry.
- ``shank`` (optional for single-shank probes): dimension: (1 x ``nSites``).
  Shank ID for each site.
  For example, ``shank = [1, 1, 1, 1, 2, 2, 2, 2]`` will assign sites 1-4 to shank 1 and
  sites 5-8 to shank 2.
  This corresponds to the :ref:`shankMap` parameter.
- ``maxSite`` (optional): scalar.
  Number of neighboring sites in each direction over which to extract spike waveforms.
  This corresponds to the :ref:`nSiteDir` parameter.
- ``nSites_ref`` (optional): scalar.
  Total number of reference sites to exclude for feature extraction.
  This corresponds to the :ref:`nSitesExcl` parameter.
- ``um_per_pix`` (optional): scalar.
  Micrometers per center-to-center vertical site spacing.
  This corresponds to the :ref:`umPerPix` parameter.

.. _io-raw:

Raw recording files
~~~~~~~~~~~~~~~~~~~

JRCLUST primarily supports the `SpikeGLX <https://github.com/billkarsh/SpikeGLX>`_ recording format,
namely a flat binary file containing 16-bit signed integers with a corresponding
`metadata <https://github.com/billkarsh/SpikeGLX/blob/master/Markdown/Metadata.md>`_ file.
ADC samples from channels 1 to :math:`n` sampled at time :math:`k` are stored together in series.
For example, let :math:`s_{i,j}` denote the sample from the :math:`i` th channel at time step :math:`j`.
Then the samples are ordered like so:

:math:`s_{1,1}`, :math:`s_{2,1}`, ..., :math:`s_{n,1}`, :math:`s_{1,2}`, :math:`s_{2,2}`, ..., :math:`s_{n,2}`, ... :math:`s_{n,k}`

For more information, see the `SpikeGLX documention on GitHub <https://github.com/billkarsh/SpikeGLX/blob/master/Markdown/UserManual.md>`__.

JRCLUST will ask for a path to your meta file when :ref:`bootstrapping <bootstrapping>` and
use this to infer a path to your raw recording.
If you don't have a meta file, don't panic!
You will be able to select your raw recording, at which point you may specify the relevant information before you begin sorting.

.. _io-trial:

Trial file
~~~~~~~~~~

This can be either a MAT-file or a CSV.
If a MAT-file, it should contain an array called ``times``.
If a CSV file, it should contain a single column.
In either case, the array or the column should contain trial start times (in seconds).
These can be used to plot the peristimulus time histogram (PSTH).

.. Additionally, you may generate this file from one of the recorded analog channel containing TTL pulses.
.. To generate the trial file from a recorded channel, run ``jrc maketrial myparam.prm``
.. You will need to supply the channel number (starting with 1), and whether to detect rising or falling edge.
.. If you have multiple recordings, the trial timing from all your files will be concatenated.

Output files
------------

+-------------------------+-------------------------------------------+---------------------------------+
| Filename                | Content                                   | Format                          |
+=========================+===========================================+=================================+
| mysession\_res.mat      | :ref:`io-res`                             | MATLAB data                     |
+-------------------------+-------------------------------------------+---------------------------------+
| mysession\_filt.jrc     | :ref:`io-filt-traces`                     | Binary (16-bit                  |
|                         | from a subset of channels                 | signed integer)                 |
+-------------------------+-------------------------------------------+---------------------------------+
| mysession\_raw.jrc      | :ref:`io-raw-traces`                      | Binary (16-bit                  |
|                         | from a subset of channels                 | signed integer)                 |
+-------------------------+-------------------------------------------+---------------------------------+
| mysession\_features.jrc | :ref:`io-features`                        | Binary (single-precision float) |
+-------------------------+-------------------------------------------+---------------------------------+
| mysession.csv           | :ref:`io-cluster-data`                    | CSV                             |
+-------------------------+-------------------------------------------+---------------------------------+

.. _io-res:

Detection and clustering results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following data are stored in mysession\_res.mat.

+-----------------------+------------------------------+------------------------+
| Variable name         | Content                      | Data format            |
+=======================+==============================+========================+
| spikeTimes            | Spike timing in ADC          | ``nSpikes`` x 1: int32 |
|                       | sample unit                  |                        |
+-----------------------+------------------------------+------------------------+
| spikeSites            | Site with the peak           | ``nSpikes`` x 1: int32 |
|                       | spike amplitude              |                        |
|                       | (center site)                |                        |
+-----------------------+------------------------------+------------------------+
| spikeSites2           | Site with the second         | ``nSpikes`` x 1: int32 |
|                       | peak spike amplitude         |                        |
+-----------------------+------------------------------+------------------------+
| spikeAmps             | Spike amplitude              | ``nSpikes`` x 1: int16 |
|                       | (local min. after            |                        |
|                       | filtering)                   |                        |
+-----------------------+------------------------------+------------------------+
| spikesBySite          | Cell of the spike            | Cell of vector of      |
|                       | indices per site             | int32                  |
+-----------------------+------------------------------+------------------------+
| siteThresh            | Detection threshold          | 1 x ``nSpikes``: single|
|                       | per site                     |                        |
+-----------------------+------------------------------+------------------------+
| filtShape             | Dimensions of                | 1 x 3: double          |
|                       | :ref:`io-filt-traces`        |                        |
+-----------------------+------------------------------+------------------------+
| rawShape              | Dimensions of                | 1 x 3: double          |
|                       | :ref:`io-raw-traces`         |                        |
+-----------------------+------------------------------+------------------------+
| featuresShape         | Dimensions of                | 1 x 3: double          |
|                       | :ref:`io-features`           |                        |
+-----------------------+------------------------------+------------------------+
| hClust                | Clustering data              | DensityPeakClustering  |
|                       |                              | object                 |
+-----------------------+------------------------------+------------------------+
| hCfg                  | Configuration data           | Config                 |
|                       |                              | object                 |
+-----------------------+------------------------------+------------------------+



.. _io-filt-traces:

Filtered spike traces
~~~~~~~~~~~~~~~~~~~~~

A binary file containing filtered spike traces, extracted in a spatiotemporal window around the spiking event.
These are stored as an
:math:`n_{\text{samples}} \times n_{\text{sites}} \times n_{\text{spikes}}`
array of signed (16-bit) integers.

:math:`n_{\text{samples}}` is determined from :ref:`evtWindow` and
:math:`n_{\text{sites}}` is determined  from :ref:`nSiteDir`.

.. _io-raw-traces:

Raw spike traces
~~~~~~~~~~~~~~~~

A binary file containing raw spike traces, extracted in a spatiotemporal window around the spiking event.
These are stored as an
:math:`n_{\text{samples}} \times n_{\text{sites}} \times n_{\text{spikes}}`
array of signed (16-bit) integers.

:math:`n_{\text{samples}}` is determined from :ref:`evtWindowRaw` and
:math:`n_{\text{sites}}` is determined  from :ref:`nSiteDir`.

.. _io-features:

Spike features
~~~~~~~~~~~~~~

A binary file containing spike features.
For each spike, JRCLUST computes some number of features at one or more positions in a
spatiotemporal window.
Consequently, the spike features data is stored as an
:math:`n_{\text{features}} \times n_{\text{positions}} \times n_{\text{spikes}}`
array of single-precision (32-bit) floating point values.

.. _io-cluster-data:

Cluster data
~~~~~~~~~~~~

A CSV export of spike-cluster data.
For each spike, a row is written containing:
- spike time
- cluster ID
- center site ID
