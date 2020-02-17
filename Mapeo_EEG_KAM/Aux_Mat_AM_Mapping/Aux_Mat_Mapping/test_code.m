% clear all;
close all;
cd('E:\Willeke2\simulation');
three_sig_simulation;
figure,plot(sig');
%% parameters
p=5;
UC=10^-3;
ksmoother=100;
fs=256;
freq=0:0.1:30;
%%
cd('E:\Willeke2\CPU_code');
[iADTF,ffADTF,swADTF,spectrum_tot]=PvM_data_to_iADTF_ffADTF_swADTF(sig',p,UC,ksmoother,fs,freq);
figure,imagesc(iADTF);
%%
cd('E:\Willeke2\GPU_code');
nr_chan=size(sig,1);
nr_samples=size(sig,2);
error=PvM_Kalman_GPU_save_memory(sig',p,UC,ksmoother,'test');
filename='test_p5_UC0.001_ks100';
[spectrum_tot]=PvM_TVAR_to_WATF_par(nr_chan,nr_samples,p,fs,freq,filename);
file=fopen([filename '_WATF_par_freq_0_to_30']);
data=fread(file,'float32');
swADTF_GPU=reshape(data,nr_chan*nr_chan,length(data)/(nr_chan*nr_chan))';