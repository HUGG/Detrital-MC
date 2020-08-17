Input data file formats
=======================

Here we explain the formats of data files that can be read for Detrital MC.

In addition to the input file (``input/det_mc_input.txt``), Detrtial MC can read observed/measured age data and three different formats of predicted age data.
Each are described in more detail below. 

Please check the rest of the documentation for more detailed explanations of how other parts the software operate.

Observed age data file format
-----------------------------

Predicted age data file format
------------------------------

Pecube Comparison.txt file
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: none
    :caption: Format of the Pecube Comparison.txt file for Detrital MC

            NLINES
         LONGITUDE    LATITUDE        ELEV PECUBE_ELEV   PECUBE_VZ     AHE_OBS    AHE_PRED     AFT_OBS    AFT_PRED     ZHE_OBS    ZHE_PRED     ZFT_OBS    ZFT_PRED     KAR_OBS    KAR_PRED     BAR_OBS    BAR_PRED     MAR_OBS    MAR_PRED     HAR_OBS    HAR_PRED   FTL_OBS01   FTL_OBS02   FTL_OBS03   FTL_OBS04   FTL_OBS05   FTL_OBS06   FTL_OBS07   FTL_OBS08   FTL_OBS09   FTL_OBS10   FTL_OBS11   FTL_OBS12   FTL_OBS13   FTL_OBS14   FTL_OBS15   FTL_OBS16   FTL_OBS17  FTL_PRED01  FTL_PRED02  FTL_PRED03  FTL_PRED04  FTL_PRED05  FTL_PRED06  FTL_PRED07  FTL_PRED08  FTL_PRED09  FTL_PRED10  FTL_PRED11  FTL_PRED12  FTL_PRED13  FTL_PRED14  FTL_PRED15  FTL_PRED16  FTL_PRED17   RAMAN_OBS  RAMAN_PRED
    ...

CSV based on Pecube Comparison.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lat,lon,elev,Pecube elev,edot,AHe,ZHe,AFT,ZFT,KAr,BAr,MAr,HAr,Mean FTL,Raman,GeolID,GlacID,MorID,RGID,Ksn,Ksn_t045,Ksn_t2,Ksn_t3,Ksn_t2m,Ksn_t3m,Ksn_t2t4,Ksn_t3t4,Ksn_t2mt4,Ksn_t3mt4,ssp,ssp_t2b31,ssp_t3b42,slope,slope-gt-30-deg,slope-lt-10-deg


Newer CSV based on Pecube Comparison.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~