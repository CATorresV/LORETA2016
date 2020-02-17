clear all;
close all;
fs=256;
t=-5:1/fs:(5-1/fs); %5s
t_seiz=0:1/fs:(5-1/fs);
t0=find(t==0);
%nr_sig=3;
SNR_brain=5; % SNR of 1/f-noise vs sinusoid
delay=[5 8];
f_seiz=[12 8];
f_sinus=linspace(f_seiz(1),f_seiz(2),length(t_seiz));

sinus=sin(2*pi*(f_sinus.*t_seiz));
sig_noise(1,:)=PvM_SimulateEEG(fs,t(1),t(end)); %1/f-noise
sig_noise(2,:)=PvM_SimulateEEG(fs,t(1),t(end));
sig_noise(3,:)=PvM_SimulateEEG(fs,t(1),t(end));
std_noise=mean(std(sig_noise'));

std_sig_wanted=std_noise*10^(SNR_brain/20);
sig=sig_noise;
sig(1,t0:end)=sinus/std(sinus)*std_sig_wanted+sig(1,t0:end);
sig(2,t0+delay(1):end)=[zeros(1,delay(1)) sig(1,t0+delay(1):end-delay(1))]+sig(2,t0+delay(1):end);
sig(3,t0+delay(2):end)=[zeros(1,delay(2)) sig(1,t0+delay(2):end-delay(2))]+sig(3,t0+delay(2):end);

%%
y=sig'; %data
p=10; % model order
UC=10^-3; %update coefficient
kalman_smoother=100; 
freq=0:0.1:30; 
[iADTF,ffADTF,swADTF,spectrum_tot]=PvM_data_to_iADTF_ffADTF_swADTF(y,p,UC,kalman_smoother,fs,freq);
