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
