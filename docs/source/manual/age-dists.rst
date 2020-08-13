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

In order to calculate age distributions for all grain ages in a sample, the first step is to calculate the probability distribution function :math:`\mathrm{PDF}(x)` for a single age assuming a normal distribution of error about the mean age :math:`\mu` with the standard deviation :math:`\sigma`, and with a kernel width scaling factor :math:`\alpha`.
This PDF calculation applies to both measured and predicted thermochronometer ages. though the predicted PDFs can be scaled by various factors, as described below.

.. math::

   \mathrm{PDF}(x) = \frac{1}{\alpha \sigma_{i} \sqrt{2 \pi}} \exp \left(-\frac{1}{2} \left(\frac{x - \mu}{\alpha \sigma} \right)^{2} \right)

Scaling for predicted age PDFs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The predicted age PDFs are scaled by one or more scaling factors :math:`f_{\mathrm{eff}}` in order to account for factors that might increase the probability of an age being present in a catchment predicted age distribution, such as differences in the tectonic uplift rate or bedrock mineral fertility.
Thus, the age PDF for a given predicted age can be calculated as

.. math::

   \mathrm{PDF}_{\mathrm{p}}(x) = f_{\mathrm{eff}} \times \mathrm{PDF}(x)

The values for the scaling factors that are combined as :math:`f_{\mathrm{eff}}` are given in the input data file for Detrital MC.
Further detail about this is given in the :doc:`section describing the Detrital MC input file <input-file>`.

Catchment sample SPDFs
----------------------

Catchment cumulative distributions
----------------------------------

Smoothed distributions
~~~~~~~~~~~~~~~~~~~~~~

Unsmoothed distributions
~~~~~~~~~~~~~~~~~~~~~~~~
