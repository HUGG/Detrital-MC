Using Detrital MC with Pecube
=============================

Introduction
------------

This is a brief overview of the steps involved in predicting detrital thermochronometer age distributions using Pecube. In general, the steps are to first predict ages across the surface, extract them for the different catchments under consideration, and then calculate predicted age distributions using different scaling factors for age prevalence in the catchments such as the instantaneous exhumation rate from Pecube, various topometric values or the abundance of ages in the different bedrock units in the catchment. The steps are detailed further below.

Step 1. Forward model cooling ages in Pecube
--------------------------------------------

**Input(s)**: Best-fit Pecube model parameters from bedrock data inversion

**Output(s)**: Predicted thermochronometer ages across the entire model surface for the different model geometries

- Surface age distributions are calculated in Pecube using the best-fit input values from inversion of the bedrock thermochronometer ages in Bhutan. Various model geometries were used depending mainly on the extent of the drainage area for the different catchments.

Step 2. Generation of catchment topometric and geological/glacial datasets
--------------------------------------------------------------------------

**Input(s)**: SRTM DEMs, processed TRMM rainfall data, bedrock geology, glacial coverage, moriane locations

**Output(s)**: Calculated normalized channel steepness index ($`k_{\mathrm{sn}}`$) and specific stream power (ssp) grid files; grids for the bedrock geology and current glacial, moraine and rock glacier coverage in the catchment

- Bodo uses the SRTM digital elevation data and his code(s) to calculate the "topometric" values for each catchment. These calculations are mainly the channel steepness index ($`k_{\mathrm{sn}}`$) and specific stream power (ssp) values. Different reference values for the reference concavity index ($`\theta`$) are used for the channel steepeness index value calculations, which is why there are several options for that. The different values for the channel steepness index and specific stream power are listed in detail under **Step 4**, where the datafiles for detrital age calculation are produced.
- The bedrock geology data is from a simplified geological map processed by Bodo to label the different bedrock units with ID numbers. The glacial coverage is from existing ArcGIS data and moraines were selected by hand using Google Earth imagery and processed by Bodo.

Step 3 - Merging topometric data
--------------------------------

**Input(s)**: Topometric, geology and glacial coverage datafiles

**Output(s)**: Merged datafiles for each catchment containing all of the data above, if applicable

- The topometric and bedrock geology/glacial files are separate for each catchment. There is one file containing all of the topometric data, a file for the bedrock geology and separate files for the glaciers, moraines and rock glaciers in each catchment, if they exist.
- These datafiles needed to be merged and this was done using a Python script written by Dave (`merge_catchment_data.py`).
    - The data point locations in each catchment were not identical, so they had to be shifted to be able to merge the files.
    - Once shifted, common points for the bedrock geology and topometric datafiles were written to a new file.
    - If there were glaciers, moraines or rock glacier in the catchment, these data were also written to the new "merged" file.

Step 4 - Adding ages from Pecube to the merged file
---------------------------------------------------

**Input(s)**: Predicted ages at surface from Pecube, merged topometric/geology/glacial datafiles

**Output(s)**: Ages and topometric datafiles for each catchment for input to **Detrital MC**

- Step 4 involves combining the "merged" topometric and geology/glacial datafiles with the predicted thermochronometer ages from Pecube. The Python code (`Pecube2catchments.py`) for doing this does the following:
    - Reads in the predicted age and instantaneous exhumation rate data from the appropriate Pecube age output file (`Ages###.vtk`) for the selected set of catchments
    - Creates an interpolation function for the Pecube elevations, instantaneous exhumation rates and predicted ages
    - Reads in the "merged" topometric/geology/glacial data file
    - Uses the interpolation functions created above to go point-by-point through the locations of data in the “merged” data file and probe for the Pecube elevation, instantaneous exhumation rate and thermochronometer age(s) at each point
    - Writes a line to the new output file containing the predicted age, exhumation rate, topometric, geology and glacial data for each point in the catchment
    - The first line of the new output file lists the number of points in the catchment. All remaining lines give the following data, separated by commas:
        - Column 1: Latitude ($`^\{\circ}`$E)
        - Column 2: Longitude ($`^\{\circ}`$N)
        - Column 3: SRTM elevation (m)
        - Column 4: Pecube interpolated elevation (m)
        - Column 5: Pecube interpolated instantaneous exhumation rate (mm/a)
        - Column 6: Pecube interpolated AHe age (Ma)
        - Column 7: Pecube interpolated ZHe age (Ma)
        - Column 8: Pecube interpolated AFT age (Ma)
        - Column 9: Pecube interpolated ZFT age (Ma)
        - Column 10: Pecube interpolated K-Ar age (Ma)
        - Column 11: Pecube interpolated B-Ar age (Ma)
        - Column 12: Pecube interpolated M-Ar age (Ma)
        - Column 13: Pecube interpolated H-Ar age (Ma)
        - Column 14: Pecube interpolated mean fission-track length ($`\mu`$m?)
        - Column 15: Pecube interpolated Raman peak T ($`^\{\circ}`$C)
        - Column 16: Bedrock geology ID
        - Column 17: Glacier ID
        - Column 18: Moraine ID
        - Column 19: Rock glacier ID
        - Column 20: Ksn, steepness index calculated with catchment-adjusted concavity values
        - Column 21: Ksn, steepness index calculated with a concavity (theta) of 0.45 (standard value)
        - Column 22: Ksn, weighted with TRMM2B31 (theta adjusted for each catchment)
        - Column 23: Ksn, weighted with TRMM3B42 (theta adjusted for each catchment) 
        - Column 24: Ksn, weighted with mean TRMM2B31 (compare this to Column 22) 
        - Column 25: Ksn, weighted with mean TRMM3B42 (compare this to Column 23) 
        - Column 26: Ksn, weighted with TRMM2B31 and theta of 0.45
        - Column 27: Ksn, weighted with TRMM3B42 and theta of 0.45
        - Column 28: Ksn, weighted with mean TRMM2B31 and theta of 0.45 (compare to Column 26)
        - Column 29: Ksn, weighted with mean TRMM3B42 and theta of 0.45 (compare to Column 27)
        - Column 30: ssp, specific stream power
        - Column 31: ssp, specific stream power weighted with TRMM2B31
        - Column 32: ssp, specific stream power weighted with TRMM3B42
- These data files are used as input predicted age files for **Detrital MC**

Step 5. - Calculating detrital age distributions with Detrital MC
-----------------------------------------------------------------

**Input(s)**: Ages and topometric data files from Step 4

**Output(s)**: Predicted and observed catchment age PDF and CDF/ECDF data files for plotting; percentage of Monte Carlo predicted age distributions that pass the Kuiper test, meaning they are statiscically indistinguishable from the observed age distribution at the $`\alpha = 0.05`$ significance level. In other words, there is a 95% probability that the two distributions are the same.

- The output files from Step 4 are part of the inputs for this step, where **Detrital MC** is used to calculate predicted and observed age distributions for each catchment. Age distribution outputs include predicted and observed normalized probability distribution functions (PDFs), cumulative distribution functions (CDFs) or empirical cumulative distribution functions (ECDFs).
- The general methodology of the operations of **Detrital MC** follows that described in [*Ruhl and Hodges* (2005)](https://doi.org/10.1029/2004TC001712) and [*Stock et al.* (2006)](https://doi.org/10.1130/G22592.1). Those steps are described in some detail below.
    - **Detrital MC** first reads the input file `input/det_mc_input.txt`. The input file defines overall what **Detrital MC** will do and the basic operating parameters. Documentation of those parameters is given in the input file.
    - **Detrital MC** then starts a loop over all catchments listed in the input file
    - In most cases, the user has selected to output the observed age PDF or CDF/ECDF, or requested comparison of the predicted and observed ages, so **Detrital MC** will:
        - Read in the observed ages for the catchment from the observed age file in `data/observed_ages/`.
        - Store the mean/median/standard deviation in the observed age uncertainties
        - Read in the predicted ages and age abundance scalings
            - Common options here inlcude scaling age prevalence in the catchment age distribution by the instantaneous exhumation rate in Pecube, bedrock outcrop in the catchment, glacial coverare or various topometric values
            - The percentage uncertainty from the observed age data will be applied to the predicted ages, unless the user provides their own percentage uncertainty
        - Create the observed age PDF and/or CDF/ECDF
            - To create the observed age PDF, PDFs are first created for all individual observed ages assuming a normal distribution of error
            - The individual age PDFs are then summed and normalized so that the area beneath the summed PDF curve is equal to 1
                - *Ruhl and Hodges* (2005) refer to this as the synoptic probability density function (SPDF)
                - This is the output PDF for the observed ages
            - If CDF or ECDF output is requested, then
                - A cumulative distribution (CDF) is produced by integrating the area beneath the observed age PDF for the catchment. This is the same data portrayed another way, and better for the statistical testing done later.
                - If ECDF output is requested, the uncertainties will not be considered in the output cumulative distribution, which is an alternative approach advocated by [*Vermeesch* (2007)](https://doi.org/10.1029/2006JF000671).
        - Unless the user has requested output of the PDF/CDF/ECDF for the entire set of catchment predicted ages, the next operation in **Detrital MC** is to start the loop for Monte Carlo sampling of the predicted basin age distribution and comparison to the data
        - The Monte Carlo loop does the following:
            - Sets the number of ages *n* to be used for calculation of the given predicted age PDF/CDF/ECDF, typically the number of ages in the observed age dataset
            - For each iteration, the code will
                - Randomly select $`n`$ ages from the predicted age distribution for the entire catchment
                - Create a predicted age PDF and/or CDF/ECDF for the selected ages following the approach described above for the observed ages
                - Compare the predicted CDF or ECDF to the data using Kuiper's test and record whether or not the two distributions were equal
        - After running over all iterations, typically 10,000, the code will write out the selected number of PDFs and/or CDFs/ECDFs in the desired format
        - The summary of the statistical comparison, most notably the percentage of Monte Carlo predicted age distributions that were statistically equal to the observed age distribution at the $`\alpha = 0.05`$ significance level will be written to the file called `kuiper_mc_results_###_samples_BASIN-NAME.dat'`
        - Finally, the code will write out the observed age PDF and/or CDF/ECDF to a file
    - If there are additional basins to process, the process above will be repeated until the end of the list of basins is reached

Notes
-----

1. Bodo notes the following about the topometric data:
	- Column 21 is likely the one you want to use if you want to compare your results to other studies
	- Column 20 is the more appropriate way of deriving steepness values
