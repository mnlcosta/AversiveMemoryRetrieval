This repository contains the scripts to reproduce the main and supplementary results of the paper: 

## Contents
•	run01_Behaviour: Reproduces the behavioral results shown in Fig. 1.
•	run02_TimeFreqAmygdalaStat: Reproduces the time-frequency results and statistics for patients with electrodes implanted in the amygdala (Fig. 2). Also produces results for Supplementary Fig. 4 (n = 17).
•	run02_TimeFreqHippocampusStat: Reproduces the time-frequency results and statistics for patients with electrodes implanted in the hippocampus (Fig. 2 and Supplementary Fig. 3, n = 12).
•	run03_GlobalERSAmygdala: Reproduces the global encoding-retrieval similarity (ERS) analysis in the amygdala (n = 17).
•	run03_GlobalERSAmygdala_control: Reproduces the control global ERS analysis in the amygdala using different frequency bands (n = 17).
•	run03_GlobalERSHippocampus: Reproduces the global ERS analysis in the hippocampus (n = 12).
•	run03_GlobalERSHippocampus_control: Reproduces the control global ERS analysis in the hippocampus using different frequency bands (n = 12).
•	run03_GlobalERSCortex: Reproduces the global ERS analysis in the lateral temporal cortex (n = 12).
•	run03_GlobalERS_amyhctcortxcomparison: Compares global ERS results in the amygdala, hippocampus, and lateral temporal cortex (n = 12).
•	run04_EESERSamygdalagammapeaks: Computes encoding-retrieval similarity using amygdala encoding patterns around amygdala gamma peaks, and tests for reactivation in the hippocampus at encoding and retrieval.
•	run04_EESERS_amyamy: Computes encoding-retrieval similarity using amygdala encoding patterns around amygdala gamma peaks, and tests for reactivation in the amygdala at retrieval.
•	run04_EESERSamygdalagammarndnopeak: Performs a control analysis using random, non-peak selections, as reported in Supplementary Fig. 12.
•	run04_EESERShippogammapeaks: Computes encoding-retrieval similarity using hippocampal encoding patterns around hippocampal gamma peaks, and tests for reactivation in the hippocampus at retrieval.
•	run05_HippocampalERSusingmaxAHEES: Computes hippocampal gamma activity patterns at encoding that are most similar to amygdala emotional encoding activity, and tests for reactivation in the hippocampus during emotional retrieval.

## Data Access
Due to size constraints, the data are not hosted in this repository. You can access all data files through the following UPM Drive folder: https://drive.upm.es/s/RygiPJJwaYBCAlf

## Requirements
- MATLAB R2019b or later
- FieldTrip toolbox
- Custom functions included in ‘/utils’
