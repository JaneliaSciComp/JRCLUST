# JRCLUST

JRCLUST is a scalable and customizable package for spike sorting on [high-density silicon probes](https://www.nature.com/articles/nature24636).
It is written in MATLAB and CUDA.

JRCLUST was originally developed by [James Jun](https://www.simonsfoundation.org/team/james-jun/) and is currently maintained by [Vidrio Technologies](https://vidriotechnologies.com).

## Installing JRCLUST

If you'd like to
test the latest development code, you can [clone the
repository](https://help.github.com/articles/cloning-a-repository/) to
your computer. If you want to stay on a release, head to the [releases
page](https://github.com/JaneliaSciComp/JRCLUST/releases) and download
the latest release.

Run the following command in MATLAB (you may want to add it to your [startup
script](https://www.mathworks.com/help/matlab/ref/startup.html)):

```matlab
addpath('/path/to/JRCLUST');
```

You may also need to recompile your CUDA codes if you're not on Windows.
Do this with

```matlab
jrclust.CUDA.compileCUDA();
```

## Questions?

Read the documentation [here](https://jrclust.readthedocs.io/en/latest/index.html) and the original bioRxiv paper [here](https://www.biorxiv.org/content/early/2017/01/30/101030).
Still have questions?
Visit our [Gitter community](https://gitter.im/JRCLUST/community) and ask!
