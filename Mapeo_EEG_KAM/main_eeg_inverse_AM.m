% Kernel-based mapping beamformer vs loretta
%__________________________________________________________________________
% Copyright (C) 2015 
% Andres Marino Alvarez Meza
% $Id: main_eeg_inverse_AM.m 
clear all; close all; clc;
%% Main toolboxes
addpath(genpath('Aux_Mat_AM_Mapping'))
addpath(genpath('FastEMD'))
%% Read transformation matrix model:
load L %lead field matrix -> conductivity matrix
M = full(L); clear L
load QG %spatial relationships
load head_model %brain model isotropic
%% Emulating a EEG
inter = [0 1];%time interval
Fs = 250; %sample frequency
t = linspace(inter(1),inter(2),Fs);
SNR = 5;%signal to noise ratio
Ndip = 1;%# of simulated sources
dipamp = [1000; 1000; 1000; 500].*1e-9; %dip amplitudes
dipfreq = [15 20 10 25];%dip frequencies	
%simulating synthetic EEG -> random localization
fprintf('Simulating EEG...')
[y,x] = eeg_stationary_simulation(t,SNR,Ndip,M,QG,vert,face,dipamp,dipfreq);
fprintf('done\n')

%% LORETA-based source localization
fprintf('Loretta-based solution...')
inverse1.type = 'LOR'; %Loreta
inverse1.type_cov = 'lineal'; %Loreta
[J0,CAY] = spm_eeg_invert_redAM(y,M,QG,vert,face,inverse1);
Y0=M*J0;%EEG reconstruction from loreta-based solution
fprintf('done\n')
%% Kernel-based solution
fprintf('kernel-based solution...')
inverse.type = 'LOR';%'LOR';%'EBB'; %Loreta
inverse.type_cov = 'corr';%'rbf'; %Loreta
%inverse.type_cov_alpha = 2;
[Jk,CAYk] = spm_eeg_invert_redAM(y,M,QG,vert,face,inverse);
Yk=M*Jk;%EEG reconstruction from loreta+kernel-based solution
fprintf('done\n')
%% plot EEG
close all
fprintf('Plotting EEG...')
figure
subplot(2,2,1)
plot(t,x,'b')
legend({'Clean EEG'})
subplot(2,2,3)
plot(t,y,'r')
legend({['Noisy EEG SNR=' num2str(SNR) '[dB]']})
subplot(2,2,2)
plot(t,Y0,'r')
legend({'Loretta solution'})
subplot(2,2,4)
plot(t,Yk,'r')
legend({'Kernel solution'})



fprintf('done\n')

%plot covariance matrices
figure
imagesc(CAY), colorbar
title('lineal')
figure
imagesc(CAYk), colorbar
title(inverse.type_cov)

%plot source energies
D = full(QG);
D = D-min(min(D));
D = D./max(max(D));
D = 1-D;
extra_mass_penalty= -1;

figure
subplot(311)
ex = sum(x.^2,2);
ex = ex/max(ex);
stem(ex);
title('Original')

subplot(312)
e0 = sum(J0.^2,2);
e0 = e0/max(e0);
stem(e0);
title(['LOR - ' inverse1.type_cov])

subplot(313)
ek = sum(Jk.^2,2);
ek = ek/max(ek);
stem(ek)
title([inverse.type '-' inverse.type_cov])

%plot source energies on the brain model
Plot3DBrain(ex,vert,face,QG); colorbar
title('Clean source')

Plot3DBrain(e0,vert,face,QG); colorbar
reo = norm(ex-e0);
title(['LOR - ' inverse1.type_cov 'MSE=' num2str(reo)])

Plot3DBrain(ek,vert,face,QG); colorbar
rek = norm(ex-ek);
title([inverse.type '-' inverse.type_cov 'MSE=' num2str(rek)])

showfigs_c(3)
%% 
%save LORLORK1
%load LORLORK1