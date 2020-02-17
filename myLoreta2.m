function [X_s] = myLoreta2( y_obs, x_prior, L,M)
% Computing the inverse problem solution following the Tikhonov
% Regularization method with LORETA for estimating the solution.
% ||Ax-y||^2 + ||Gamma x||^2
% y_obs size N x M
% X_s   size P x M
% A     size N x P
n=size(x_prior,1);
d=size(y_obs,1);
%Gamma in Tikhonov regularization is defined as Gamma = alpha * I
alpha = 0.4; %Regularization parameter
lambda = alpha * eye(n);
% 
% X_s = inv(A' * A + Gamma' * Gamma) * A' * y_obs ;

% dir ='C:\Users\Alejandro\Desktop\Problema inverso\SPM\spm12\Datasets\';
% file= 'sim_2016cdbespm12_SPM_CTF_MEG_example_faces1_3D.mat';
% load([dir file]);
% A=struct(D);
% clear D
% se = A.sensors.meg.chanpos;
% clear A

%% Laplacian matrix operator construction


% L= zeros(n,n);
% 
% df=4; % reference distance
% 
% for i=1:n
%     for j=1:n
%       
%         if (i==j)
%             L(i,j)=1;
%         elseif (sqrt(i^2-j^2)<=df)
%             L(i,j)=-1/6;
%         end
%     end
% end


% A_0 = (1/2)*( eye(M) + inv( eye(M) * diag(A_1 * ones(M))))*A_1;
% 
% A = kron(A_0,eye(3));
% 
% B = (6/d^2) * (A - eye(3*M));
%% Non- modeled features

Sigma = cov(y_obs');

Ce = Sigma^2 * eye(d); %Covariance matrix


%% Lead Field Matrix
% a= linspace(-100,100,floor(nthroot(n,3)));
% b= linspace(-100,100,floor(nthroot(n,3))+1);
% c= linspace(-60,60,floor(nthroot(n,3))+1);
% [x, y, z] = ndgrid(a,c,b);
% 
% V = [x(:),y(:),z(:)];
% V = V(1:n,:);
%V=LF;
%N = 132; %number of electrodes
%% Electrodes positions
%S=se;%str2num(importdata('positions.txt'));
% [x,y,z]=sphere;
% x=(x(:));
% y=(y(:));
% z=(z(:)); %Centered at (0,0,0)
% 
% x=x((z>0));
% y=y((z>0));
% z=z((z>0));
% 
% x=x(1:d);
% y=y(1:d);
% z=z(1:d);
% 
% S = [x,y,z];   %Coordinates of the measurement electrodes

%K = zeros(N,M,3); %Initialization of the lead field matrix
%% MAtrix of relationships

% sigma = 0.3; %Conductivity of the medium
% M = zeros(d,n);
% 
% for i=1:d
%     for j=1:n
%         term1 = (1/(4*pi*sigma))*((S(i,:)-V(j,:))/(norm(S(i,:)-V(j,:),3)));
%         term2 = (1/(4*pi*sigma))*((S(i,:)-V(j,:))/(norm(S(1,:)-V(1,:),3)));
%         k=mean(term1-term2);
%         M(i,j) = k;
%     end
% end

clear V
clear LF


%% LORETA
fprintf('computing A...');
A= (M' /Ce) * M + lambda^2 * (L' * L);
fprintf('computing B...');
B=((M' / Ce) *y_obs + lambda^2 * (L' * L) * x_prior);
fprintf('Solving Loreta ... ');
X_s = A\B;

fprintf('... Done \n')
%% PAra visualización

% X_c= (X_s-min(X_s))/(max(X_s)-min(X_s)) ;
% 
% X_c = im2uint8(X_c);
% X_c(X_c==0)=1;
% 
% r = linspace(1,255,256);
% g = ones(1,256) * 100;
% b = linspace(-255,1,256) * -1;
% color = zeros(length(X_c),3);
% 
% for i=1:length(X_c)
%      color(i,:)= [round(r(X_c(i))) round(g(X_c(i))) round(b(X_c(i)))];
% end
% 
% color=color/255;
% figure
% scatter3(V(:,1),V(:,2),V(:,3),20,color,'fill') 
% colormap summer
% hold on
% scatter3(S(:,1),S(:,2),S(:,3),40,'r')
