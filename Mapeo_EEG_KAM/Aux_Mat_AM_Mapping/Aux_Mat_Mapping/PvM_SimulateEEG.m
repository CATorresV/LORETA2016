function eeg = PvM_SimulateEEG(fs,ti,tf)
% fs = sampling frequency
% ti = start time
% tf = end time
% This function generates a simulated background eeg with some superimposed
% brain activity in the alfa band (8-13 Hz)


%% Prameters initialization
% Define the parameters necessary for the simulation
%clear; clc
noe = 1; %number of epochs
noc = 1; %number of channels
%fs = 256; %sampling frequency
T = 1/fs; %sampling interval
%ti = 0; % starting time in sec
%tf = 1; % ending time in sec
t = ti:T:tf; %vettore asse t per un'epoca
%alfa = 1;
%display(sprintf('PARAMETERS \n'));
%display(sprintf('T = %2f sec', T));
%display(sprintf('fs = %2d Hz', fs));
%display(sprintf('ti = %2d sec', ti));
%display(sprintf('tf = %2d sec', tf));
%display(sprintf('alfa = %2d \n', alfa));

%% Creation of the random input with a flat spectrum
% First a vector of random numbers normally distributed, according to
% N(0,1) is generated, then the obtained signal is filtered in the
% frequency band (0.1 Hz-100 Hz). The simulated scalp eeg is supposed to
% have an amplitude between 10 microV and 50 microV (this means zero mean
% and 25 std
%display(sprintf('CREATION OF THE RANDOM INPUT \n'));
%fprintf('Creating and filtering the random input...');
x = 25*randn(noc,length(t)*noe);
xfilt = filtro(x,0.1,99,fs,20);
%display(sprintf('Done \n'));

%% Calculate the spectrum of the random input

%fprintf('Computing the fourier transform...');
X = fft(xfilt);
N = length(X);
%calcolo delle frequenze
v = 1/N/T; %int di freq 
if mod(N,2)==0
    f = (-N/2+(0:N-1))*v;
else f = (-(N-1)/2+(0:N-1))*v;
end;
%display(sprintf('Done \n'));


%% Modification of the spectrum
% We modified the spectrum to obtain a 1/f behaviour: to do that we
% multiplied the original spectrum by a function OOF = 1/f

%fprintf('Modifying the spectrum...');
halfOOF = 1./(f(find(f==0):end)+1);

if mod(N,2)==0
    OOF = [halfOOF fliplr(halfOOF)];
else OOF = [halfOOF fliplr(halfOOF(2:end))];
end

Xmod = OOF.*X;
%display(sprintf('Done \n'));


%% Calculate the signal in the time domain
% After the modification of the spectrum, we transformed the signal back to
% the time domain to obtain the simulated background eeg and we computed the
% power spectrum to check the presence of the 1/f behaviour: we expected a
% slope of 2 in the PSD plot (due to the square of the 1/f function)

%fprintf('Calculating the time domain signal...');

xmod = ifft(Xmod,'symmetric');
%the signal is rescalming to obtain a physiologic amplitude
xmod = (xmod-mean(xmod))*std(xfilt)/std(xmod);
[Pxfilt,ffilt] = pwelch(xfilt,[],[],128,fs,'onesided');
[Pxmod,fmod] = pwelch(xmod,[],[],128,fs,'onesided');
eeg=xmod;
%display(sprintf('Done \n'));

%display(sprintf('mean xmod = %2f ', mean(xmod)));
%display(sprintf('std xmod= %2f ', std(xmod)));


