.. _migrating:

Migrating from JRCLUST v3
=========================

What's different
----------------

Variable and parameter names have changed.
The input and output file names have changed.
A probe file is no longer required, but may be used.
vcFile and csFile_merge are now a single parameter, rawRecordings.
Some commands are missing.
Logging and curation history are stored in the ``hClust`` object instead of in a separate file.
You must specify your parameter file for every invocation of JRCLUST.
