Spike detection
===============

Summary
-------
.. include:: ../detect-summary.rst

.. _denoising:

Denoising
---------

JRCLUST constructs an adaptive notch filter which cleans narrow-band noise peaks exceeding a threshold you supply :ref:`threshold you supply <fftThreshMAD>`.
(Set this parameter (``fftThreshMAD``) to 0 to skip this step.)

The time-domain signal is converted to the frequency domain via fast Fourier transform (FFT).
After this, the average power across channels is computed.
Power is (putatively) detrended using :math:`\frac 1 f` relationship in frequency vs. power.
For each frequency bin, the power-frequency product is normalized by its `MAD`_, and those frequencies whose MAD power-frequency products exceed this threshold are zeroed out.
The cleaned frequency-domain signal is then transformed back to the time domain via inverse FFT.

The figure below shows the original (left) and detrended (right) plots.
Red indicates FFT coefficients that exceed threshold or points that are within 3 samples from the outlier points.

.. image:: /.static/fftClean.png

.. _filtering:

Filtering
---------

This step filters the raw samples depending on the :ref:`filter type <filterType>` you specify.
Additionally, JRCLUST computes the common average reference and subtracts it from the filtered samples.
(The meaning of "average" in "common average reference" is controlled by the :ref:`CAR mode <carMode>` you specify.)

.. _compute-threshold:

Computing detection thresholds
------------------------------

.. A threshold is automatically computed from the filtered samples in the following manner due to
.. `Quiroga <https://www.mitpressjournals.org/doi/abs/10.1162/089976604774201631>`__
.. :

Assuming a Gaussian noise distribution on every site, we estimate the standard
deviation :math:`\sigma_{\text{noise}}^{(i)}` of the background on site :math:`i` by

.. math::

    \DeclareMathOperator{\med}{med}
    \sigma_{\text{noise}}^{(i)} = \med\left( \frac{|x^{(i)}|}{0.6745} \right)

where :math:`x^{(i)}` is the filtered signal on site :math:`i`.
(See `here <http://www.scholarpedia.org/article/Spike_sorting#Step_ii.29_Spike_Detection>`__ for justification.)
A single :ref:`multiple <qqFactor>` of :math:`\sigma_{\text{noise}}^{(i)}` that you specify
serves as the detection threshold on each site.
In other words, each site has a distinct threshold which is a multiple of an estimate of
the standard deviation of its own noise distribution.

.. _peak-detection:

Detecting peaks
---------------

Samples on each site are compared to the threshold :math:`\text{Thr}_i` for that site.
Those samples with negative value but exceeding the threshold in magnitude are considered
putative peaks.
(If you :ref:`specify <fDetectBipolar>`, then JRCLUST can also detect peaks with positive value.)
A putative peak must also exceed its immediate neighbors in magnitude, i.e., the sample at time
:math:`t_j` is a peak if and only if the following conditions hold:

.. math::

    \DeclareMathOperator{\myAmp}{amp}
    \begin{align*}
        |\myAmp(t_j)| &> \text{Thr} \\
        |\myAmp(t_j)| &> |\myAmp(t_{j-1})| \\
        |\myAmp(t_j)| &> |\myAmp(t_{j+1})|
    \end{align*}

You may also impose :ref:`additional constraints <nneigh_min_detect>` on peak detection.

.. _merge-peaks:

Merging peaks
-------------

Peak samples are compared with neighbors in a spatiotemporal window around each peak to detect duplicates.
If peaks are found within this window, then the following operations are performed:

#. If a neighboring peak has a larger amplitude, then discard the center peak.
#. If a neighboring peak has the same amplitude but occurs first, then discard the center peak.
#. If a neighboring peak has the same amplitude and occurs at the same time, then discard the peak
   occurring on the higher-indexed site.

.. _extract-windows:

Extract spatiotemporal windows around spiking events
----------------------------------------------------



.. _MAD: https://en.wikipedia.org/wiki/Median_absolute_deviation
