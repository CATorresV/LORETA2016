% clear all;
close all;
fs=256;
t=-5:1/fs:(5-1/fs); %5s
t_seiz=0:1/fs:(5-1/fs);
t0=find(t==0);
nr_sig=2;
SNR_brain=5; % SNR of 1/f-noise vs sinusoid
delay=5;
f_seiz=[12 8];
f_sinus=linspace(f_seiz(1),f_seiz(2),length(t_seiz));

sinus=sin(2*pi*(f_sinus.*t_seiz));
sig_noise(1,:)=PvM_SimulateEEG(fs,t(1),t(end));
sig_noise(2,:)=PvM_SimulateEEG(fs,t(1),t(end));
std_noise=mean(std(sig_noise'));

std_sig_wanted=std_noise*10^(SNR_brain/20);
sig=sig_noise;
sig(1,t0:end)=sinus/std(sinus)*std_sig_wanted+sig(1,t0:end);
sig(2,t0+delay:end)=[zeros(1,delay) sig(1,t0+delay:end-delay)]+sig(2,t0+delay:end);
