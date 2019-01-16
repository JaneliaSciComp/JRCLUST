Spike detection
===============

   The sample at time :math:`t` is a peak if and only if the following conditions hold:

.. math::

    \DeclareMathOperator{\myAmp}{amp}
    \begin{align*}
        |\myAmp(t)| &> thresh \\
        |\myAmp(t)| &> |\myAmp(t-1)| \\
        |\myAmp(t)| &> |\myAmp(t+1)|
    \end{align*}
