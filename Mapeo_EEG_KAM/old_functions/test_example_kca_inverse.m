%test_example
clear; close; clc;

load L
M = full(L); clear L
load QG
load head_model

inter = [0 1];
Fs = 250;
t = linspace(inter(1),inter(2),Fs);

SNR = 7;


Ndip = 3;
dipamp = [1000; 1000; 1000].*1e-9;
dipfreq = [15 20 10];	

[y,x] = eeg_stationary_simulation(t,SNR,Ndip,M,QG,vert,face,dipamp,dipfreq);


%Initial solution (LORETA)
inverse.type = 'LOR'; %Loreta
J0 = spm_eeg_invert_red(y,M,QG,vert,face,inverse);
Y0=M*J0;

sp = kScaleOptimization_info_alpha(squareform(pdist(y)),2);
% sp = 1;
D = squareform(pdist(y));
Ky = exp(-D.^2/(2*sp^2));

etav = [1e-3 1e-3];
J  = kMetric_Learning_Mahalanobis_2A_s(M,Ky,ones(size(M,1),1),size(y,2),true,etav,J0,sp);
Ykca = M*J;
close all

%EMD
D = full(QG);
D = D-min(min(D));
D = D./max(max(D));
D = 1-D;
extra_mass_penalty= -1;

figure
subplot(311)
stem(sum(x.^2,2));
subplot(312)
stem(sum(J0.^2,2));
subplot(313)
stem(sum(J.^2,2))

% EMDy(1)=emd_hat_gd_metric_mex(mean(x.^2,2),mean(J0.^2,2),D,extra_mass_penalty);
% EMDy(2)=emd_hat_gd_metric_mex(mean(x.^2,2),mean(J.^2,2),D,extra_mass_penalty);
% 
Plot3DBrain(sum(x.^2,2),vert,face,QG);
Plot3DBrain(sum(J0.^2,2),vert,face,QG);
% title(num2str(EMDy(1)))
Plot3DBrain(sum(J.^2,2),vert,face,QG);
% title(num2str(EMDy(2)))


