
function [yn,x]=eeg_stationary_simulationA(t,SNR,Ndip,M,QG,vert,face,dipamp,dipfreq,randiA)
%forward eeg simulation
%yn:noisy output
%x:clean states
%tint:time interval
%Eduardo Giraldo
%egiraldos@utp.edu.co

%carga matrices de campo guia y de vecinos
%simulation Juan David

Nd = size(M,2);
Nc = size(M,1);

startf1  = 0;					% (sec) start time
duration = 0.6;					% duration 

endf1 = duration + startf1;
f1ind = intersect(find(t>startf1),find(t<=endf1));

% -----------------

signal = zeros(Ndip,length(t));
for j=1:Ndip				% For each source
		f1 = dipfreq(j);	% Frequency depends on stim condition
		amp1 = dipamp(j);	% also the amplitude
		phase1 = pi/2;
		signal(j,f1ind) = signal(j,f1ind)...
			+ amp1*sin((t(f1ind)...
			- t(min(f1ind)))*f1*2*pi + phase1);
end

act_dip		= zeros(Ndip,1);
act_dipo		= zeros(Ndip,3);
act_dip(1)		= randiA;%randi(Nd);
act_dipo(1,:)	= vert(act_dip(1),:);
if Ndip>2
    for i=2:Ndip
        act_dip(i)		= randi(Nd);
        act_dipo(i,:)	= vert(act_dip(i),:);
        dists			= zeros(i-1,1);
        for j = 1:i-1
            dists(j)	= pdist2(act_dipo(i,:),act_dipo(j,:));
        end
        while sum(dists<10) ~= 0
            act_dip(i)		= randi(Nd);
            act_dipo(i,:)	= vert(act_dip(i),:);
            for j = 1:i-1
                dists(j)	= pdist2(act_dipo(i,:),act_dipo(j,:));
            end
        end
    end
end

X	  = zeros(Nd,length(t));						% Matrix of dipoles
% Add waveform of all smoothed sources to their equivalent dipoles
% QGs add up to 0.9854
dists = dist([vert; act_dipo]' );
act_dip = [];
for i = 1:Ndip
    [~, act_dip(i)] = min(dists(end-i+1,1:end-Ndip));
end
for j=1:Ndip
    for i=1:length(t)
        X(:,i) = X(:,i) + signal(j,i)*QG(:,act_dip(j));
    end
end
Y = M*X;    % Forward model

allchanstd=(std(Y'));
meanrmssignal=mean(allchanstd);
MEANCHANNOISE=1; %% TAKE MEAN RMS NOISE RATHER THAN A DIFFERENT NOISE LEVEL PER CHANNEL
for i = 1:size(Y,1)
    if Ndip>0,
        if ~MEANCHANNOISE,
            channoise = std(Y(i,:)).*randn(size(Y(i,:)))/(10^(SNR/20));
        else
            channoise = meanrmssignal.*randn(size(Y(i,:)))/(10^(SNR/20));
        end;
    else
        channoise = randn(size(Y(i,:))); %% just add noise when no dips specifed
    end
    allchannoise(i,:)=channoise;
    Y(i,:) = Y(i,:) + channoise;
end
yn = Y;
x = X;
