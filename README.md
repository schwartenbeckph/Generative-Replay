# Generative-Replay

Code to run the main analyses from https://www.biorxiv.org/content/10.1101/2021.06.06.447249v1.

Note that here there's only individual preprocessed MEG data for the first two subjects due to space limitations. Data from all subjects can be found at https://www.dropbox.com/sh/c52di7x1tt2k3jp/AABPov3W6-c9dH0TViBECHj7a?dl=0.

Replay results can be reproduced based on summary stats in results file 'Replay_InferenceTime.mat' and 'Replay_InferenceLag.mat' or re-run using individual pre-processed data. Either way, start with the 'Replay_InferenceTime.m' (Figures 6 D-F) and 'Replay_InferenceLag.m' (Figure 6C) script.

MEG RSA can be reproduced using 'RSA_inference.m' (Figure 5). Summary stats 'RSA_inference.mat' can be found in link above.

Second level univariate fMRI results can be found in 'fMRI_SummaryStats_Univariate' in the above link. Code used to produce the univariate results (Figure 2, Figure 3C) is in 'Stats_FirstLevel_Suppression.m' in the fMRI_code folder.

Second level multivariate ('silhouette algebra') fMRI results can be found in 'fMRI_SummaryStats_Algebra' in the above link. Code used to produce the silhouette algebra results (Figure 3B) is in 'mk_RSA_STIMS_ALGEBRA_ols_vol.m' in the fMRI_code folder.
