# Supplementary

This document provides information about the data that was used for the paper titled “Accounting for Spatial Anonymization in DHS Household Surveys.

1.) Survey Data

This paper is constructed mainly on the data set of Kenya 2014 Standard Demographic Health Survey. The survey is conducted by The Demographic and Health Surveys (DHS) Program (https://dhsprogram.com/data/dataset_admin/index.cfm).

Access to the data set requires registration. This can be done by following the steps that are explained in the  web page of DHS (https://dhsprogram.com/data/Using-DataSets-for-Analysis.cfm#CP_JUMP_14037).

*Following files are used for analyses:

-Survey locations are obtained from a shape file called KEGE71FL.shp. 
-Survey responses are obtained from a file called KEIR72DT/KEIR72FL.DTA. Answers to the question "Ever used anything or tried to delay or avoid getting pregnant”are available as the variable v302a of this data set. Types of the contraceptive methods  related to these answers are available as variable v312. Detailed explanation of the variables in the data sets are available in the document called “Standard Recode Manual for DHS-7”. (ICF. 2018. Demographic and Health Surveys Standard Recode Manual for DHS7. The Demographic and Health Surveys Program. Rockville, Maryland, U.S.A.: ICF)


2.) Administrative Area (Level-1) Borders 

Level-1 administrative area borders are needed to see the county that each cluster center is located before and after jittering. If the borders will be respected while jittering, the code will test the newly created location against the true one and keep jittering the true location until both of them are within the same county. 

This data file can be obtained from https://gadm.org/index.html.

If there is a need for converting the coordinates from degrees to kilometers or vice versa, there are functions for both operations in functions.R script.
