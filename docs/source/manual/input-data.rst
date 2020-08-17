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
               LON         LAT        ELEV PECUBE_ELEV   PECUBE_VZ     AHE_OBS    AHE_PRED     AFT_OBS    AFT_PRED     ZHE_OBS    ZHE_PRED     ZFT_OBS    ZFT_PRED     KAR_OBS    KAR_PRED     BAR_OBS    BAR_PRED     MAR_OBS    MAR_PRED     HAR_OBS    HAR_PRED   FTL_OBS01   FTL_OBS02   FTL_OBS03   FTL_OBS04   FTL_OBS05   FTL_OBS06   FTL_OBS07   FTL_OBS08   FTL_OBS09   FTL_OBS10   FTL_OBS11   FTL_OBS12   FTL_OBS13   FTL_OBS14   FTL_OBS15   FTL_OBS16   FTL_OBS17  FTL_PRED01  FTL_PRED02  FTL_PRED03  FTL_PRED04  FTL_PRED05  FTL_PRED06  FTL_PRED07  FTL_PRED08  FTL_PRED09  FTL_PRED10  FTL_PRED11  FTL_PRED12  FTL_PRED13  FTL_PRED14  FTL_PRED15  FTL_PRED16  FTL_PRED17   RAMAN_OBS  RAMAN_PRED
               ...

CSV file format 1
~~~~~~~~~~~~~~~~~

.. code-block:: none
    :caption: Format of CSV file 1 for Detrital MC

    NLINES
    UNUSED,UNUSED,UNUSED,UNUSED,SCALE1,AGE_PRED1,AGE_PRED2,AGE_PRED3,AGE_PRED4,UNUSED,UNUSED,AGE_PRED5,UNUSED,UNUSED,UNUSED,UNIT_ID1,UNIT_ID2,UNIT_ID3,UNIT_ID4,SCALE2,SCALE3,SCALE4,SCALE5,UNUSED,UNUSED,UNUSED,UNUSED,UNUSED,UNUSED,UNUSED,SCALE6,SCALE7


CSV file format 2
~~~~~~~~~~~~~~~~~

.. code-block:: none
    :caption: Format of CSV file 2 for Detrital MC

    NLINES
    UNUSED,UNUSED,UNUSED,UNUSED,SCALE1,AGE_PRED1,AGE_PRED2,AGE_PRED3,AGE_PRED4,UNUSED,UNUSED,AGE_PRED5,UNUSED,UNUSED,UNUSED,UNIT_ID1,UNIT_ID2,UNIT_ID3,UNIT_ID4,SCALE2,SCALE3,SCALE4,SCALE5,UNUSED,UNUSED,UNUSED,UNUSED,UNUSED,UNUSED,UNUSED,SCALE6,SCALE7,SCALE8,UNIT_ID5,UNIT_ID6

