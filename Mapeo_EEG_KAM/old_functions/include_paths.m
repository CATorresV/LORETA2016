%Paths
clear; close; clc;

%spm 12
addpath('D:\KCA inverse\spm12b\spm12b');
spm('defaults','EEG');

%FastEMD - Error measure
addpath(genpath('FastEMD'));

%Tools
addpath('simulation')

%KCA
addpath('D:\KCA inverse\CODES KCA\MLmat-master_Relevancia');