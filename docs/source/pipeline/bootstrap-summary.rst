Your config file encapsulates the choices in parameters that you make, as well as describing the relevant probe geometry.
You can create a config file by invoking ``jrc bootstrap`` from the MATLAB command window.

This will collect recording-specific information from your `SpikeGLX meta file (.meta) <https://github.com/billkarsh/SpikeGLX/blob/master/Markdown/Metadata.md>`__ as well as setting default parameters.
The default parameters may be inspected on the `parameters`_ page.

If you don't have a .meta file or your .meta file is missing some data, JRCLUST requests the following information:

* **Sample rate (Hz)**: read from ``imSampRate`` (or ``niSampRate`` for NI recordings) in .meta file
* **\# channels saved**: read from ``nSavedChans`` in .meta file
* **Microvolts per bit (uV/bit)**: computed from ``imAiRangeMax``, ``imAiRangeMin`` (or ``niAiRangeMax``, ``niAiRangeMin``, and ``niMNGain`` for NI recordings) in .meta file
* **Header offset (bytes)**: set to 0 for SpikeGLX recordings since no header is stored in the .bin file
* **Data type**: select from ``int16``, ``uint16``, ``single``, or ``double`` (SpikeGLX files are saved as ``int16``)
* **Neuropixels probe option**: 0 if N/A; read from ``imroTbl`` in .meta file

.. _`parameters`: ../parameters/index.html
