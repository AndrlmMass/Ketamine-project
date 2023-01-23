t = addpath('Z:/Data/Anesthesia_Project/Kongsberg_ED_Study/psych/');
addpath('Z:/Scripts/templates/eeg_channel_templates/BrainVisionN/active electrodes/_ARCHIV/actiCAP/actiCAP 64 Channel/')
savepath('Z:/Data/Anesthesia_Project/Kongsberg_ED_Study/psych/pathdef.m')

EEG = pop_loadbv(t, 'SD5001_AWAKE_REST_EC.vhdr');

% Import the BFEV file
EEG = eeg_checkset(EEG, 'load', 'AC-64.bfev');

% Change the channel locations using the BFEV file as the template
EEG = pop_chanedit(EEG, 'changefield',{1 'X'}, 'AC-64.bfev');

