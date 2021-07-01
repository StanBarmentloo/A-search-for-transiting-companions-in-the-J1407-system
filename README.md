# J1407_transit_search_activity
Repository for the code used to create the figures in the paper 'A search for transiting companions in the J1407 system'.

Timeseries data for the six telescopes considered in the paper is stored in the 'Data' directory. 

Arguably the most important data product produced by our research is the combined, 20 year base line ground based telescope activity removed dataset. This data is created in the 'Figure_2' directory and stored as 'Final_Combined_Data.csv' in the 'Data' directory. The data is sorted per MJD, but the 'Telescope' column allows to separate it per telescope. The numbers in this column correspond to the following telescopes: 1 = ASAS, 2 = ASAS-SN, 3 = KELT, 4 = PROMPT, 5 = ROAD. Similar results for the TESS data can be found as 'Final_TESS_Data.csv' in the 'Data' directory and are created in the 'Figure_3' directory. 

The 'Figure_x' directories are self consistent (together with the data from 'Data') and reproduce each figure from the paper via Jupyter notebooks and Python scripts.
