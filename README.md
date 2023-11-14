# finetrack_prep
Routine to read an ATL03 file and output calibration data and expected waveform tables

Output h5 file (atl03_calibrations_and_wf_tables.h5) is to be used as input when computing fine_track heights using python fine_track routine, which requires expected waveform distributions and CAL19 data.

Data contained within atl03_calibrations_and_wf_tables.h5:

-atl03_filename
-orientation
-CAL19/
--dead_time
--fpb_corr
--strength
--width
-exp_wf_table/
--binz
--mz
--sdz
--wf_table1
--wf_table2
-gtx/
--dead_time
--strong_weak
--wf_table_select


Sample execution:

python finetrack_prep.py -g1 ATL03_20220726163210_05311604_006_02.h5


