SEIZURE DURATION AND FREQUENCY ...
Code and data for figure generation and statistical analysis

The code provided on this repository generates results and supplementary figures for - publication, as well as the corresponding statistical analysis.

# CODE
The provided code, gmb_SNS_ASM_Dur_V7.py, has been developed with Spyder v6 running on python 3.12.5

- Structure
This code is structured for modular execution, allowing the generation of each figure independently by running specific sections. NOTE: the way section are assigned might not be interpreted the same on other platforms.

The code can be executed all at once to generate and store all figures in both svg and pdf format.

In case of being interested on a specific figure, please follow this steps. 
1. Execute the section titled "Initialise Workspace", to generate and load all required data. 
2. Execute the section corresponding to the figure of interest, regardless of order on the file. 

- Configuration
The data used to generate the results has to be placed on a subfolder named "data" in the same folder as the main code is stored.
The variable 'sz_Clust'  is the main source of information for the presented results and supplementary analysis. It is loaded from the following css file: '/data/rednmf_l10.csv'. Check line 36 to modify the origin of the data.
The variable 'Stt_Freq_asm' contains the seizure frequency information. It loaded from the following css file: '/data/rednmf_l10.csv'. Check line 37 to modify the origin of the data.

The code is configured to automatically generate figures and store them as both .svg and .pdf figures. 
To disable figure storage, change 'save_fig' variable (line 40) to "False"
'path_save_fig' variable identifies the root file to store figures when 'save_fig' has been set to "True". 
The selected interpreter, 'cairo', does not show figures on Sypder interface. To use standard figure interpreter, please comment lines 25 to 30.

#DATA
Two files are provided as the main data source for figure generation. All the provided data has been extracted considering ROI segmentation on Lausanne scale 60, and band power on base 10 logarithmic scale. For different configuration, please contact the authors.

- rednmf_l10_szFrq.csv
This file contain seizure frequency for all individual included on the full dataset (28 individuals). Seizure frequency has been calculated on the two ASM periods defined on this study: Normal-dose and Under-dose.

- rednmf_l10.csv
This file contains the main data for figure generation and statistical analysis. Each row corresponds to a seizure and each column to a variable. The included variables relate to: Patient identification, seizure states, asm-related seizure states reassignment, duration (full, under-dose state duration and common sates), number of states (full, under-dose and common), seizure classification, ASM levels and state, time of day (including circular 24h and 12h varialbles), type of ASM in use for tapering. 

As explained on the publication, the initial complete dataset, composed of 28 patients and +300 seizures, is reduced as more restrictive analysis is applied (from 28 patients to 17 and finally 7). As a result several variables present empty, nan and ±inf values. This data has been taken into account an removed before any figure generation or statistical analysis.
