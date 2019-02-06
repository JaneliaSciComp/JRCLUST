.. _migrating:

Migrating from JRCLUST v3
=========================

What's different
----------------

- Variable and parameter names have changed.
- The input and output file names have changed.
- A probe file is no longer required, but may be used.
- ``vcFile`` and ``csFile_merge`` are now a single parameter, ``rawRecordings``.
- Some commands are :ref:`missing <missing-commands>`.
- Logging and curation history are stored in the DensityPeakClustering object ``hClust``, instead of in a separate file where they might get out of sync.
- You must specify your parameter file for every invocation of JRCLUST.
- JRCLUST will overwrite old config files and save a backup copy to 'filename.prm.bak'.

.. note::

	Old sortings must be manually imported with ``jrc importv3 mysession_jrc.mat``.
	This will output ':ref:`mysession_res.mat <io-res>`' and update your config file.
	(A backup will be made of your config file at 'mysession.prm.bak'.)
	Additionally, your '_spkraw.jrc', '_spkwav.jrc', and '_spkfet.jrc' will be renamed
	to ':ref:`_raw.jrc <io-raw-traces>`', ':ref:`_filt.jrc <io-filt-traces>`', and
	':ref:`_features.jrc <io-features>`', respectively.
	Nothing will be altered except for this renaming.
