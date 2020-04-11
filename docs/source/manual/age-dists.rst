Calculation of age distributions
================================

This page presents and overview of how age distributions in Detrital MC are calculated. Distributions of ages from detrital samples can be assembled and visualized in serveral different ways. Below I describe the different distributions and their meanings, as well as how they are used in Detrital MC.

Clearly, the page is under construction.

.. admonition:: Terminology and abbreviations

   - **Probability density function (PDF)**: 
   - **Synoptic probability density function (SPDF)**:
   - **Cumulative density function (CDF)**:
   - **Empirical cumulative density function (ECDF)**:

Individual age PDFs
-------------------

In order to calculate age distributions for all grain ages in a sample, the first step is to calculate the probability distribution function :math:`\mathrm{PDF}_{i}` for a single age assuming a normal distribution of error about the mean age :math:`\bar{x}_{i}` with the standard deviation :math:`\sigma_{i}`.

.. math::

   \mathrm{PDF}_{i} = \frac{1}{\sigma_{i} \sqrt{2 \pi}} \exp \left(-\frac{1}{2} \left(\frac{x - \bar{x}_{i}}{\sigma_{i}} \right)^{2} \right)

Catchment sample SPDFs
----------------------

Catchment cumulative distributions
----------------------------------

Smoothed distributions
~~~~~~~~~~~~~~~~~~~~~~

Unsmoothed distributions
~~~~~~~~~~~~~~~~~~~~~~~~
