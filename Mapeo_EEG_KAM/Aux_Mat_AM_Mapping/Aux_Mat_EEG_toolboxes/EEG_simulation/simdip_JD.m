%simdip_JD
%
% Created by:	Jose David Lopez - ralph82co@gmail.com
%
% Date: January 2012
% Set of possible dipoles: 1930, 3714, 5859
% Dipole for bad location: 710

function D = simdip_JD(D)

if isfield(D.inv{1},'inverse')
	D.inv{1} = rmfield(D.inv{1},'inverse');
end

%% GENERATION OF NEURAL SOURCES
% SPM works in mm, lead fields are correctly calculated in SI
% units (spm_cond_units)
% Select one, two, or three sources

% One source
meshsourceind = 3714;				% Surce position locator in ctf space
Ndip = size(meshsourceind,2);		% Number of dipoles
dipfreq = 20;						% Source frequency
dipamp = 1000.*1e-9;				% Source amplitude in nAm

% Two synchronous sources
% meshsourceind = [5859 7044];		% Surce position locator in ctf space
% Ndip = size(meshsourceind,2);		% Number of dipoles
% dipfreq = [15 4];					% Source frequency
% dipamp = [1000; 1000].*1e-9;		% Source amplitude in nAm

% Three sources
%[7044 417 4419];%
% meshsourceind = [5859 1930 3714];	% Surce position locator in ctf space
% Ndip = size(meshsourceind,2);		% Number of dipoles
% dipfreq = [15 20 10];				% Source frequency
% dipamp = [1000; 1000; 1000].*1e-9;	% Source amplitude in nAm

%% WAVEFORM FOR EACH SOURCE

Ntrials = D.ntrials;				% Number of trials

% define period over which dipoles are active
startf1  = 0;					% (sec) start time
duration = 0.6;					% duration 

endf1 = duration + startf1;
f1ind = intersect(find(D.time>startf1),find(D.time<=endf1));

% Create the waveform for each source
signal = zeros(Ndip,length(D.time));
for j=1:Ndip				% For each source
	for i=1:Ntrials			% and trial
		f1 = dipfreq(j);	% Frequency depends on stim condition
		amp1 = dipamp(j);	% also the amplitude
		phase1 = pi/2;
		signal(j,f1ind) = signal(j,f1ind)...
			+ amp1*sin((D.time(f1ind)...
			- D.time(min(f1ind)))*f1*2*pi + phase1);
	end
end


%% CREATE A NEW FORWARD PROBLEM

fprintf('Computing Gain Matrix: ')
[L D] = spm_eeg_lgainmat(D);				% Gain matrix
Nd    = size(L,2);							% number of dipoles
X	  = zeros(Nd,161);						% Matrix of dipoles
fprintf(' - done\n')

% Green function for smoothing sources with the same distribution than SPM8
fprintf('Computing Green function from graph Laplacian:')
vert  = D.inv{1}.mesh.tess_mni.vert;
face  = D.inv{1}.mesh.tess_mni.face;
A     = spm_mesh_distmtx(struct('vertices',vert,'faces',face),0);
GL    = A - spdiags(sum(A,2),0,Nd,Nd);
GL    = GL*0.6/2;
Qi    = speye(Nd,Nd);
QG    = sparse(Nd,Nd);
for i = 1:8
    QG = QG + Qi;
    Qi = Qi*GL/i;
end
QG    = QG.*(QG > exp(-8));
QG    = QG*QG;
clear Qi A GL
fprintf(' - done\n')

% Add waveform of all smoothed sources to their equivalent dipoles
% QGs add up to 0.9854
for j=1:Ndip 
	for i=1:161
		X(:,i) = X(:,i) + signal(j,i)*QG(:,meshsourceind(j));
	end
end

D(:,:,1) = L*X;					% Forward problem

% Copy data to all trials (It should not be the same, of course)
for i=1:Ntrials
	D(:,:,i) = D(:,:,1);
end

%% ADD WHITE NOISE
% In decibels

SNR = 10;		% Select SNR in dB
for i = 1:size(D,1)
	channoise = std(D(i,:,1)).*randn(size(D(i,:,1)))/(10^(SNR/10));
	D(i,:,1) = D(i,:,1) + channoise;
end

%% Plot and save

graf_SPM_JD(D,X);

figure
hold on
aux = L(1,:)*X;
plot(D.time,aux,'r');
hold on
plot(D.time,D(1,:,1));
title('Measured activity over MLC11');
legend('Noiseless','Noisy');

D.save;

fprintf('\n Finish\n')

